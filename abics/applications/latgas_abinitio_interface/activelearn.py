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

"""
Building blocks of the reference-data step of active learning (``abics_mlref``).

The functions here operate on one configuration or one file at a time.
They do not use MPI, do not change the working directory and do not exit
the process; every directory is passed explicitly.  ``abics_mlref``
(:mod:`abics.scripts.activelearn`) composes them for the staged workflow,
and an integrated driver can call them directly.

Directory layout (``rootdir`` is where ``abics_mlref`` is run)::

    rootdir/ALloop.progress                 stage sequence (AL0, train, MC0, ...)
    rootdir/rundirs.txt                     solver input directories to be run
    rootdir/AL<n>/baseinput.progress        finished relaxation stages of AL<n>
    rootdir/AL<n>/<replica>/input<i>/initial.vasp
                                            structure before perturbation
    rootdir/AL<n>/<replica>/input<i>/baseinput<stage>/
                                            solver input and output
    rootdir/AL<n>/<replica>/structure.<i>.vasp
                                            final structure returned by the solver
    rootdir/AL<n>/<replica>/energy_corr.dat MLIP energy, reference energy, index
    rootdir/AL<n>/<replica>/relax_max.dat   largest displacement per configuration
    rootdir/MC<n>/<replica>/                output of abics_sampling
"""

from __future__ import annotations

import glob
import os
from typing import Iterable, Iterator, Optional, Sequence, Tuple

from pymatgen.core import Structure

from abics.applications.latgas_abinitio_interface.model_setup import (
    perturb_structure,
)

PROGRESS_FILE = "ALloop.progress"
STAGE_FILE = "baseinput.progress"
RUNDIRS_FILE = "rundirs.txt"


# ---------------------------------------------------------------------------
# Paths and bookkeeping files
# ---------------------------------------------------------------------------


def input_dir_format(nmax: int) -> str:
    """Format string for the names of the input directories.

    ``input_dir_format(10).format(3)`` is ``"input03"``: the index is
    zero-padded to the number of digits of ``nmax``.
    """
    ndigits = len(str(nmax))
    return "input{:0>" + str(ndigits) + "d}"


def stage_dir(input_dir: str, stage: int) -> str:
    """Solver input/output directory of relaxation stage ``stage``."""
    return os.path.join(input_dir, f"baseinput{stage}")


def list_input_indices(al_replica_dir: str) -> list[int]:
    """Indices of the ``input<i>`` directories under ``AL<n>/<replica>``."""
    sample_list = [
        int(os.path.basename(path).split("t")[1])
        for path in glob.iglob(os.path.join(al_replica_dir, "input*"))
    ]
    sample_list.sort()
    return sample_list


def read_stage(aldir: str) -> int:
    """Relaxation stage whose solver outputs are to be processed next.

    ``aldir`` is ``AL<n>``.  Returns 0 when no stage has been recorded yet.
    """
    path = os.path.join(aldir, STAGE_FILE)
    if not os.path.exists(path):
        return 0
    with open(path, "r") as fi:
        return int(fi.readlines()[-1]) + 1


def record_stage(aldir: str, stage: int) -> None:
    """Record that the solver inputs of stage ``stage`` have been processed."""
    _append_line(os.path.join(aldir, STAGE_FILE), str(stage))


def last_progress(rootdir: str) -> str:
    """Last entry of ``ALloop.progress``."""
    with open(os.path.join(rootdir, PROGRESS_FILE), "r") as fi:
        return fi.readlines()[-1].strip()


def record_progress(rootdir: str, label: str) -> None:
    """Append ``label`` (``"AL<n>"``) to ``ALloop.progress``."""
    _append_line(os.path.join(rootdir, PROGRESS_FILE), label)


def write_rundirs(rootdir: str, rundirs: Sequence[str]) -> None:
    """Write the list of solver input directories to be run to ``rundirs.txt``."""
    with open(os.path.join(rootdir, RUNDIRS_FILE), "w") as fi:
        fi.write("\n".join(rundirs))
        fi.write("\n")
        fi.flush()
        os.fsync(fi.fileno())


def _append_line(path: str, line: str) -> None:
    with open(path, "a") as fi:
        fi.write(line + "\n")
        fi.flush()
        os.fsync(fi.fileno())


# ---------------------------------------------------------------------------
# Configurations
# ---------------------------------------------------------------------------


def generate_random_structures(config, ndata: int) -> Iterator[Structure]:
    """Random configurations for the first active-learning step.

    Yields ``config.structure`` itself after each shuffle, so the caller may
    perturb it in place; the structures are generated lazily so that the
    random draws of shuffling and perturbation interleave per configuration.
    """
    for _ in range(ndata):
        constraint_fullfilled, msg = config.shuffle()
        print(msg)
        config.structure.sort(key=lambda site: site.species_string)
        yield config.structure


def count_mc_samples(mc_replica_dir: str) -> int:
    """Number of samples recorded in ``obs.dat`` of one replica."""
    nsamples = 0
    with open(os.path.join(mc_replica_dir, "obs.dat")) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            nsamples += 1
    return nsamples


def load_mc_structures(
    mc_replica_dir: str,
    sample_list: Iterable[int],
    ignore_species: Optional[Sequence[str]],
) -> Iterator[Structure]:
    """Structures sampled by ``abics_sampling``, with ignored species restored.

    ``structure.<i>.vasp`` lacks the species ignored by the sampling solver;
    they are taken from ``structure_norel.<i>.vasp`` (the structure before
    the solver run, which contains them).
    """
    for i in sample_list:
        st = Structure.from_file(os.path.join(mc_replica_dir, f"structure.{i}.vasp"))
        ignore_structure = Structure.from_file(
            os.path.join(mc_replica_dir, f"structure_norel.{i}.vasp")
        )
        if ignore_species:
            _append_ignored_species(st, ignore_structure, ignore_species)
        st.sort(key=lambda site: site.species_string)
        yield st


def _append_ignored_species(
    st: Structure, full_structure: Structure, ignore_species: Sequence[str]
) -> None:
    """Append the sites of ``ignore_species`` in ``full_structure`` to ``st``."""
    ignore_structure = full_structure.copy()
    remove_sp = filter(
        lambda sp: sp not in ignore_species, ignore_structure.symbol_set
    )
    ignore_structure.remove_species(remove_sp)
    for site in ignore_structure:
        st.append(site.species_string, site.frac_coords)


# ---------------------------------------------------------------------------
# Solver input and output of one configuration
# ---------------------------------------------------------------------------


def write_reference_input(
    solver_input, structure: Structure, perturb: float, input_dir: str, stage: int = 0
) -> str:
    """Write the solver input of the first relaxation stage for one configuration.

    ``structure`` is perturbed in place by ``perturb``; the unperturbed
    structure is saved as ``input_dir/initial.vasp``.  Returns the absolute
    path of the solver input directory.
    """
    structure0 = structure.copy()
    perturb_structure(structure, perturb)
    solver_input.update_info_by_structure(structure)
    inputfile = os.path.abspath(stage_dir(input_dir, stage))
    solver_input.write_input(inputfile)
    structure0.to(
        fmt="POSCAR",
        filename=os.path.abspath(os.path.join(input_dir, "initial.vasp")),
    )
    return inputfile


def write_next_stage_input(
    solver_input, structure: Structure, input_dir: str, stage: int
) -> str:
    """Write the solver input of relaxation stage ``stage`` from ``structure``.

    Returns the absolute path of the solver input directory.
    """
    solver_input.update_info_by_structure(structure)
    inputfile = os.path.abspath(stage_dir(input_dir, stage))
    solver_input.write_input(inputfile)
    return inputfile


def parse_reference_output(
    solver_output,
    input_dir: str,
    stage: int,
    ignore_species: Optional[Sequence[str]] = None,
) -> Tuple[float, Structure]:
    """Energy and final structure of one configuration from the solver output.

    If ``ignore_species`` is given, those species (absent from the solver
    output) are restored from ``input_dir/initial.vasp``.
    """
    energy, st = solver_output.get_results(stage_dir(input_dir, stage))
    if ignore_species:
        initial = Structure.from_file(os.path.join(input_dir, "initial.vasp"))
        _append_ignored_species(st, initial, ignore_species)
    return energy, st


def max_displacement(st_in: Structure, st_rel: Structure) -> Tuple[int, float]:
    """Largest site displacement between two structures with the same site order.

    Returns ``(site index, distance)``.
    """
    dmax = 0
    dmax_id = 0
    for j, site in enumerate(st_in):
        d = site.distance(st_rel[j])
        if d > dmax:
            dmax = d
            dmax_id = j
    return dmax_id, dmax


# ---------------------------------------------------------------------------
# Per-replica result files
# ---------------------------------------------------------------------------


def write_energy_corr(path: str, rows: Iterable[Tuple]) -> None:
    """Write ``energy_corr.dat``: MLIP energy, reference energy, sample index.

    The MLIP energy is ``"NaN"`` for the first active-learning step.
    """
    with open(path, "w") as fi:
        for ene_mlip, ene_ref, step in rows:
            fi.write(f"{ene_mlip} {ene_ref} {step}\n")


def write_relax_max(path: str, rows: Iterable[Tuple[int, float]]) -> None:
    """Write ``relax_max.dat``: site index and largest displacement per configuration."""
    with open(path, "w") as fi:
        for dmax_id, dmax in rows:
            fi.write("{} \t {}\n".format(dmax_id, dmax))
