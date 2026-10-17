# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.

"""
Check the outputs of a staged active-learning integration test.

Usage (run in the test directory after run.sh has executed
``AL.sh`` and ``MC.sh`` alternately ``nloops`` times)::

    python3 ../check_active_learn.py input.toml [nloops]

``nloops`` defaults to 2 (AL0 -> train -> MC0 -> AL1 -> train -> MC1).

The checks pin down the file-based interface between ``abics_mlref``,
``abics_train`` and ``abics_sampling`` so that a refactoring of the
active-learning scripts can be verified not to change it:

* ``ALloop.progress`` records the stage sequence.
* ``AL<n>/<replica>/`` holds ``ndata`` solver input directories,
  ``initial.vasp`` for each, the parsed structures ``structure.<i>.vasp``,
  ``energy_corr.dat`` (MLIP energy, reference energy, MC sample index) and,
  for ``n >= 1``, ``relax_max.dat``.
* For ``n >= 1`` the MLIP energies in ``energy_corr.dat`` are the ones
  recorded by ``MC<n-1>`` in ``obs.dat`` on the lines of the selected
  samples, and the reference energies are reproduced by the mock solver
  from ``structure.<i>.vasp``.
* ``MC<n>/<replica>/`` holds the sampled structures and observables.
* For ``n >= 1``, ``AL<n>/convergence.dat`` reproduces the error statistics
  (RMSE, MAE, largest error) of ``energy_corr.dat`` per replica and overall.

The reference solver of the test must be ``mock.py`` in the test directory
(``mlref.solver.function = "mock.calc_energy"``).
"""

import glob
import os
import re
import sys

import numpy as np
import toml
from pymatgen.core import Structure


class CheckFailure(Exception):
    pass


def check(cond, msg):
    if not cond:
        raise CheckFailure(msg)


def read_progress():
    with open("ALloop.progress") as f:
        return [line.strip() for line in f if line.strip()]


def count_data_lines(path):
    n = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                n += 1
    return n


def check_al_step(n, r, params, nsamples, calc_energy):
    mlref = params["mlref"]
    ndata = mlref["ndata"]
    sampler = mlref.get("sampler", "linspace")
    energy_ref = params["train"].get("energy_ref", 0.0)
    mult = int(np.prod(params["config"].get("supercell", [1, 1, 1])))

    aldir = os.path.join(f"AL{n}", str(r))
    check(os.path.isdir(aldir), f"{aldir} does not exist")

    inputdirs = sorted(glob.glob(os.path.join(aldir, "input*")))
    check(
        len(inputdirs) == ndata,
        f"{aldir}: expected {ndata} input directories, found {len(inputdirs)}",
    )
    for d in inputdirs:
        check(
            os.path.isdir(os.path.join(d, "baseinput0")),
            f"{d}/baseinput0 does not exist",
        )
        check(
            os.path.isfile(os.path.join(d, "initial.vasp")),
            f"{d}/initial.vasp does not exist",
        )

    corr_path = os.path.join(aldir, "energy_corr.dat")
    check(os.path.isfile(corr_path), f"{corr_path} does not exist")
    corr = np.loadtxt(corr_path, ndmin=2)
    check(
        corr.shape == (ndata, 3),
        f"{corr_path}: expected shape {(ndata, 3)}, found {corr.shape}",
    )
    ene_mlip = corr[:, 0]
    ene_ref = corr[:, 1]
    steps = corr[:, 2]
    check(
        np.all(steps == np.round(steps)),
        f"{corr_path}: sample indices (3rd column) must be integers",
    )
    steps = steps.astype(int)

    if n == 0:
        # No MLIP exists yet: random configurations, no prediction.
        check(
            np.all(np.isnan(ene_mlip)),
            f"{corr_path}: 1st column must be NaN in AL0",
        )
        check(
            np.array_equal(steps, np.arange(ndata)),
            f"{corr_path}: sample indices in AL0 must be 0..{ndata - 1}",
        )
    else:
        check(
            np.all(np.isfinite(ene_mlip)),
            f"{corr_path}: 1st column must be finite for AL{n}",
        )
        if sampler == "linspace":
            expected = np.linspace(0, nsamples - 1, num=ndata, dtype=np.int64)
            check(
                np.array_equal(steps, expected),
                f"{corr_path}: sample indices {steps.tolist()} differ from "
                f"linspace selection {expected.tolist()}",
            )
        else:
            check(
                len(set(steps.tolist())) == ndata
                and steps.min() >= 0
                and steps.max() < nsamples,
                f"{corr_path}: sample indices must be {ndata} distinct values "
                f"in [0, {nsamples})",
            )
        # The MLIP energy of sample i is the internal energy on line i of
        # obs.dat (the column named "energy_internal" in its header), the
        # line for which structure.<i>.vasp was written.
        obs_lines = []
        column = None
        with open(os.path.join(f"MC{n - 1}", str(r), "obs.dat")) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    name, _, label = line[1:].partition(":")
                    if label.strip() == "energy_internal":
                        column = int(name.strip().lstrip("$")) - 1
                elif line:
                    obs_lines.append(float(line.split()[column]))
        check(column is not None, f"MC{n - 1}/{r}/obs.dat: no energy_internal column in header")
        expected_mlip = np.array(obs_lines)[steps] + energy_ref * mult
        check(
            np.allclose(ene_mlip, expected_mlip, rtol=1e-10, atol=1e-8),
            f"{corr_path}: MLIP energies do not match MC{n - 1}/{r}/obs.dat",
        )

    check(
        np.all(np.isfinite(ene_ref)),
        f"{corr_path}: reference energies (2nd column) must be finite",
    )

    # abics_train adds structure.<i>_mapped.vasp next to these; count only
    # the structures written by abics_mlref.
    structures = [
        p
        for p in glob.glob(os.path.join(aldir, "structure.*.vasp"))
        if re.fullmatch(r"structure\.\d+\.vasp", os.path.basename(p))
    ]
    check(
        len(structures) == ndata,
        f"{aldir}: expected {ndata} structure.<i>.vasp files, found {len(structures)}",
    )
    for e_ref, step in zip(ene_ref, steps):
        st_path = os.path.join(aldir, f"structure.{step}.vasp")
        check(os.path.isfile(st_path), f"{st_path} does not exist")
        e_recalc = calc_energy(Structure.from_file(st_path))
        check(
            np.isclose(e_ref, e_recalc, rtol=1e-8, atol=1e-6),
            f"{corr_path}: reference energy {e_ref} for sample {step} is not "
            f"reproduced by the mock solver ({e_recalc})",
        )

    if n > 0:
        relax_path = os.path.join(aldir, "relax_max.dat")
        check(os.path.isfile(relax_path), f"{relax_path} does not exist")
        relax = np.loadtxt(relax_path, ndmin=2)
        check(
            relax.shape == (ndata, 2),
            f"{relax_path}: expected shape {(ndata, 2)}, found {relax.shape}",
        )
        # relax_max.dat is the largest displacement between initial.vasp
        # (before perturbation) and the solver's final structure.  The mock
        # solver does not relax, and perturb_structure moves every atom by
        # exactly mlref.solver.perturb, so each entry must equal that value.
        perturb = mlref["solver"].get("perturb", 0.0)
        check(
            np.allclose(relax[:, 1], perturb, rtol=0.0, atol=1e-8),
            f"{relax_path}: max displacements {relax[:, 1].tolist()} must all "
            f"equal mlref.solver.perturb = {perturb} (mock solver does not relax)",
        )


def check_convergence(n, nreplicas):
    """AL<n>/convergence.dat (n >= 1) reproduces the statistics of energy_corr.dat."""
    path = os.path.join(f"AL{n}", "convergence.dat")
    check(os.path.isfile(path), f"{path} does not exist")
    rows = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            check(len(fields) == 6, f"{path}: expected 6 columns, found {line!r}")
            rows[fields[0]] = fields[1:]
    expected_names = [str(r) for r in range(nreplicas)] + ["all"]
    check(
        list(rows) == expected_names,
        f"{path}: expected rows {expected_names}, found {list(rows)}",
    )

    all_diff = []
    for r in range(nreplicas):
        corr = np.loadtxt(os.path.join(f"AL{n}", str(r), "energy_corr.dat"), ndmin=2)
        diff = corr[:, 0] - corr[:, 1]
        all_diff.append(diff)
        _check_stats_row(path, str(r), rows[str(r)], diff, corr[:, 2].astype(int))
    all_diff = np.concatenate(all_diff)
    _check_stats_row(path, "all", rows["all"], all_diff, None)


def _check_stats_row(path, name, fields, diff, samples):
    n, rmse, mae, max_error, max_sample = fields
    check(int(n) == len(diff), f"{path}: row {name}: n = {n}, expected {len(diff)}")
    expected = {
        "rmse": np.sqrt(np.mean(diff**2)),
        "mae": np.mean(np.abs(diff)),
        "max_error": np.max(np.abs(diff)),
    }
    for key, value in zip(("rmse", "mae", "max_error"), (rmse, mae, max_error)):
        check(
            np.isclose(float(value), expected[key], rtol=1e-5, atol=1e-9),
            f"{path}: row {name}: {key} = {value}, expected {expected[key]}",
        )
    if samples is not None:
        imax = int(np.argmax(np.abs(diff)))
        check(
            int(max_sample) == samples[imax],
            f"{path}: row {name}: max_error_sample = {max_sample}, expected {samples[imax]}",
        )


def check_mc_step(n, r, nsamples):
    mcdir = os.path.join(f"MC{n}", str(r))
    check(os.path.isdir(mcdir), f"{mcdir} does not exist")
    check(
        os.path.isfile(os.path.join(f"MC{n}", "kTs.npy")),
        f"MC{n}/kTs.npy does not exist",
    )

    obs_path = os.path.join(mcdir, "obs.dat")
    check(os.path.isfile(obs_path), f"{obs_path} does not exist")
    nlines = count_data_lines(obs_path)
    check(
        nlines == nsamples,
        f"{obs_path}: expected {nsamples} samples, found {nlines}",
    )

    obs_save_path = os.path.join(mcdir, "obs_save.npy")
    check(os.path.isfile(obs_save_path), f"{obs_save_path} does not exist")
    obs = np.load(obs_save_path)
    check(
        obs.ndim == 2 and obs.shape[0] == nsamples,
        f"{obs_save_path}: expected {nsamples} rows, found shape {obs.shape}",
    )
    check(
        np.all(np.isfinite(obs[:, 0])),
        f"{obs_save_path}: energies must be finite",
    )

    for i in range(nsamples):
        for name in (f"structure.{i}.vasp", f"structure_norel.{i}.vasp"):
            check(
                os.path.isfile(os.path.join(mcdir, name)),
                f"{mcdir}/{name} does not exist",
            )


def check_rundirs(nloops, nreplicas, ndata):
    check(os.path.isfile("rundirs.txt"), "rundirs.txt does not exist")
    with open("rundirs.txt") as f:
        rundirs = [line.strip() for line in f if line.strip()]
    check(
        len(rundirs) == nreplicas * ndata,
        f"rundirs.txt: expected {nreplicas * ndata} entries, found {len(rundirs)}",
    )
    last = os.path.abspath(f"AL{nloops - 1}")
    for d in rundirs:
        check(os.path.isdir(d), f"rundirs.txt: {d} does not exist")
        check(
            os.path.abspath(d).startswith(last + os.sep),
            f"rundirs.txt: {d} is not under AL{nloops - 1}",
        )


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    params = toml.load(sys.argv[1])
    nloops = int(sys.argv[2]) if len(sys.argv) > 2 else 2

    sys.path.insert(0, os.getcwd())
    import mock  # the test's reference solver

    nreplicas = params["mlref"]["nreplicas"]
    ndata = params["mlref"]["ndata"]
    nsteps = params["sampling"]["nsteps"]
    sample_frequency = params["sampling"].get("sample_frequency", 1)
    nsamples = len(range(0, nsteps, sample_frequency))

    expected_progress = []
    for n in range(nloops):
        expected_progress += [f"AL{n}", "train", f"MC{n}"]

    try:
        progress = read_progress()
        check(
            progress == expected_progress,
            f"ALloop.progress: expected {expected_progress}, found {progress}",
        )
        for n in range(nloops):
            for r in range(nreplicas):
                check_al_step(n, r, params, nsamples, mock.calc_energy)
                check_mc_step(n, r, nsamples)
            if n > 0:
                check_convergence(n, nreplicas)
        check_rundirs(nloops, nreplicas, ndata)
    except CheckFailure as e:
        print(f"CHECK FAILED: {e}")
        sys.exit(1)

    print(
        f"CHECK OK: {nloops} active-learning loops, "
        f"{nreplicas} replicas, {ndata} reference calculations per loop"
    )


if __name__ == "__main__":
    main()
