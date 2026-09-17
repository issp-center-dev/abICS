# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2020- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

"""Generate a default set of aenet training input files (generate/train/predict)
from the [config] and [train] sections of an abICS input.toml file.

The produced files mirror the ones found in the examples
(examples/active_learning_*/aenet_train_input) and are meant as a
reasonable starting point that the user can then tune by hand.
"""

from __future__ import annotations

import argparse
import os
import sys

import toml

from abics.applications.latgas_abinitio_interface.defect import (
    DFTConfigParams,
    defect_config,
)

FINGERPRINT_TEMPLATE = """DESCR
  N. Artrith and A. Urban, Comput. Mater. Sci. 114 (2016) 135-150.
  N. Artrith, A. Urban, and G. Ceder, Phys. Rev. B 96 (2017) 014112.
END DESCR

ATOM {species}

ENV {nenv}
{env_species}

RMIN {rmin}d0

BASIS type=Chebyshev
radial_Rc = {radial_rc}  radial_N = {radial_n} angular_Rc = {angular_rc}  angular_N = {angular_n}
"""


def get_species(params_root: dict) -> list[str]:
    """Species used for training = all species in [config] minus train.ignore_species."""
    configparams = DFTConfigParams.from_dict(params_root["config"])
    config = defect_config(configparams)
    all_species = set(config.structure.symbol_set)
    ignore_species = set(params_root.get("train", {}).get("ignore_species", None) or [])
    species = sorted(all_species - ignore_species)
    if not species:
        raise ValueError(
            "No species left for training after removing ignore_species "
            f"({sorted(ignore_species)}) from config species ({sorted(all_species)})."
        )
    return species


def write_generate_input(outdir: str, species: list[str]) -> None:
    os.makedirs(outdir, exist_ok=True)
    lines = ["OUTPUT aenet.train", "", "TYPES", str(len(species))]
    for sp in species:
        lines.append(f"{sp} -0.0  ! eV")
    lines += ["", "SETUPS"]
    for sp in species:
        lines.append(f"{sp}   {sp}.fingerprint.stp")
    lines.append("")
    with open(os.path.join(outdir, "generate.in.head"), "w") as f:
        f.write("\n".join(lines) + "\n")

    for sp in species:
        env_species = "\n".join(species)
        content = FINGERPRINT_TEMPLATE.format(
            species=sp,
            nenv=len(species),
            env_species=env_species,
            rmin=RMIN,
            radial_rc=RADIAL_RC,
            radial_n=RADIAL_N,
            angular_rc=ANGULAR_RC,
            angular_n=ANGULAR_N,
        )
        with open(os.path.join(outdir, f"{sp}.fingerprint.stp"), "w") as f:
            f.write(content)


def write_train_input(outdir: str, species: list[str]) -> None:
    os.makedirs(outdir, exist_ok=True)
    lines = [
        "TRAININGSET aenet.train",
        f"TESTPERCENT {TESTPERCENT}",
        f"ITERATIONS  {ITERATIONS}",
        "",
        f"MAXENERGY {MAXENERGY}",
        "",
        "TIMING",
        "",
        "!SAVE_ENERGIES",
        "",
        "METHOD",
        "!lm batchsize=4000 learnrate=0.1d0 iter=3 conv=0.001 adjust=5.0",
        "bfgs",
        "",
        "!METHOD 2 Steepest Descent",
        "!online_sd gamma=5.0d-7 alpha=0.25d0",
        "!online_sd gamma=1.0d-8 alpha=0.25d0",
        "!",
        "!METHOD 3 Extended Kalman Filter",
        "!ekf",
        "!",
        "!METHOD 4 Levenberg-Marquardt",
        "!lm batchsize=8000 learnrate=0.1d0 iter=3 conv=0.001 adjust=5.0",
        "",
        "NETWORKS",
        "! atom   network         hidden",
        "! types  file-name       layers  nodes:activation",
    ]
    for sp in species:
        lines.append(f"  {sp}     {sp}.{HIDDEN_LAYERS_TAG}.nn    {HIDDEN_LAYERS}")
    with open(os.path.join(outdir, "train.in"), "w") as f:
        f.write("\n".join(lines) + "\n")


def write_predict_input(outdir: str, species: list[str]) -> None:
    os.makedirs(outdir, exist_ok=True)
    lines = ["TYPES", str(len(species))]
    lines += species
    lines += ["", "", "NETWORKS"]
    for sp in species:
        lines.append(f"{sp}    {sp}.{HIDDEN_LAYERS_TAG}.nn")
    lines += ["", "", "VERBOSITY low", ""]
    with open(os.path.join(outdir, "predict.in"), "w") as f:
        f.write("\n".join(lines) + "\n")


# Defaults taken from the descriptor/network settings used in the examples
# (examples/active_learning_*/aenet_train_input)
RMIN = 0.55
RADIAL_RC = 8.0
RADIAL_N = 16
ANGULAR_RC = 6.5
ANGULAR_N = 4
TESTPERCENT = 10
ITERATIONS = 500
MAXENERGY = 10000
HIDDEN_LAYERS = "2      15:tanh 15:tanh"
HIDDEN_LAYERS_TAG = "15t-15t"


def main_impl(params_root: dict, output_dir: str, force: bool) -> None:
    train_params = params_root.get("train", {})
    if train_params.get("type", "aenet") != "aenet":
        print(
            f"Warning: train.type = '{train_params.get('type')}' is not 'aenet'; "
            "generated files are for the aenet trainer only.",
            file=sys.stderr,
        )

    species = get_species(params_root)

    if output_dir is None:
        base_input_dir = train_params.get("base_input_dir", "./aenet_train_input")
        if isinstance(base_input_dir, list):
            base_input_dir = base_input_dir[0]
        output_dir = base_input_dir

    generate_dir = os.path.join(output_dir, "generate")
    train_dir = os.path.join(output_dir, "train")
    predict_dir = os.path.join(output_dir, "predict")

    for d in (generate_dir, train_dir, predict_dir):
        if os.path.exists(d) and os.listdir(d) and not force:
            print(
                f"Error: {d} already exists and is not empty. Use --force to overwrite.",
                file=sys.stderr,
            )
            sys.exit(1)

    write_generate_input(generate_dir, species)
    write_train_input(train_dir, species)
    write_predict_input(predict_dir, species)

    print(f"Generated default aenet training input files for species {species} in {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate default aenet training input files (generate/train/predict) "
            "from the [config] and [train] sections of an abICS input.toml file."
        )
    )
    parser.add_argument(
        "tomlfile", nargs="?", default="input.toml", help="abICS input file (default: input.toml)"
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory (default: train.base_input_dir from the toml file, "
        "or ./aenet_train_input)",
    )
    parser.add_argument(
        "--force", action="store_true", help="Overwrite existing non-empty output directories"
    )
    args = parser.parse_args()

    params_root = toml.load(args.tomlfile)
    main_impl(params_root, args.output_dir, args.force)


if __name__ == "__main__":
    main()
