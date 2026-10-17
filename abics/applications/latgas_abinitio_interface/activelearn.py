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
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
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


# ---------------------------------------------------------------------------
# Candidate configurations from a Monte Carlo run
# ---------------------------------------------------------------------------


@dataclass
class ReferenceCandidate:
    """A configuration sampled by ``abics_sampling`` that may be sent to the reference solver.

    This is the sampler-independent interface between the Monte Carlo output
    and the selection of reference calculations: a sampler only has to
    provide a list of these (see :func:`register_candidate_reader`).
    """

    index: int
    "Index of the sample in the replica's trajectory (used to name the input directory)"
    structure_file: str
    "Structure as seen by the sampling solver (ignored species absent)"
    structure_norel_file: str
    "Structure before the solver run, including the ignored species"
    energy: float
    "Energy of the structure predicted by the model used for sampling"
    weight: float = 1.0
    "Relative statistical weight of the sample (1 for samples of a Boltzmann distribution)"


CandidateReader = Callable[[str], List[ReferenceCandidate]]
_candidate_readers: Dict[str, CandidateReader] = {}


def register_candidate_reader(sampler_type: str, reader: CandidateReader) -> None:
    """Register how to list the candidates of a sampler.

    ``reader`` receives the replica directory (``MC<n>/<replica>``) and
    returns the candidates in trajectory order.  Samplers that use the
    output layout of ``DefaultObserver`` need no registration.
    """
    _candidate_readers[sampler_type] = reader


def reference_candidates(mc_replica_dir: str, sampler_type: str) -> List[ReferenceCandidate]:
    """Candidates sampled by ``sampler_type`` in one replica directory."""
    reader = _candidate_readers.get(sampler_type, default_candidates)
    return reader(mc_replica_dir)


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


def default_candidates(mc_replica_dir: str) -> List[ReferenceCandidate]:
    """Candidates from the output layout written by ``DefaultObserver``.

    ``obs.dat`` has one line per sample, ``structure.<i>.vasp`` and
    ``structure_norel.<i>.vasp`` are written by ``DefaultObserver.writefile``
    with the line number ``i``, and column 0 of ``obs_save.npy`` is the
    energy.  All samplers shipped with abICS (RXMC, PAMC, parallelRand,
    parallelMC) use this layout.
    """
    nsamples = count_mc_samples(mc_replica_dir)
    obs = np.load(os.path.join(mc_replica_dir, "obs_save.npy"))
    return [
        ReferenceCandidate(
            index=i,
            structure_file=os.path.join(mc_replica_dir, f"structure.{i}.vasp"),
            structure_norel_file=os.path.join(mc_replica_dir, f"structure_norel.{i}.vasp"),
            energy=obs[i, 0],
        )
        for i in range(nsamples)
    ]


def select_candidates(
    candidates: Sequence[ReferenceCandidate], refparams
) -> List[ReferenceCandidate]:
    """Choose the candidates to be calculated with the reference solver.

    ``refparams`` is the ``RefParams`` of the ``[mlref]`` section; its
    ``sampling`` method gives the positions in ``candidates`` (``ndata``
    of them, evenly spaced or random).  The candidates' weights are not
    used by these two policies.
    """
    positions = refparams.sampling(len(candidates))
    return [candidates[i] for i in positions]


def load_candidate_structures(
    candidates: Iterable[ReferenceCandidate],
    ignore_species: Optional[Sequence[str]],
) -> Iterator[Structure]:
    """Structures of the candidates, with ignored species restored.

    ``structure_file`` lacks the species ignored by the sampling solver;
    they are taken from ``structure_norel_file`` (the structure before the
    solver run, which contains them).
    """
    for candidate in candidates:
        st = Structure.from_file(candidate.structure_file)
        ignore_structure = Structure.from_file(candidate.structure_norel_file)
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


# ---------------------------------------------------------------------------
# Error statistics of the model against the reference solver
# ---------------------------------------------------------------------------


@dataclass
class ErrorStatistics:
    """Errors of the model energies against the reference energies of one AL step."""

    n: int
    "Number of configurations with both energies available"
    rmse: float
    "Root mean square error"
    mae: float
    "Mean absolute error"
    max_error: float
    "Largest absolute error"
    max_error_sample: int
    "Sample index of the configuration with the largest error (-1 if n == 0)"


def read_energy_corr(path: str) -> np.ndarray:
    """Rows of ``energy_corr.dat`` as an array of shape (n, 3)."""
    return np.loadtxt(path, ndmin=2)


def error_statistics(energy_corr: np.ndarray) -> ErrorStatistics:
    """Error statistics from rows of ``energy_corr.dat`` (model, reference, sample).

    Rows whose model energy is not available (``NaN``, as in ``AL0``) are
    skipped.  Errors are in the units of the energies, i.e. per supercell.
    """
    if len(energy_corr) == 0:
        return ErrorStatistics(0, float("nan"), float("nan"), float("nan"), -1)
    ene_model = energy_corr[:, 0]
    ene_ref = energy_corr[:, 1]
    samples = energy_corr[:, 2].astype(int)
    valid = np.isfinite(ene_model) & np.isfinite(ene_ref)
    n = int(valid.sum())
    if n == 0:
        return ErrorStatistics(0, float("nan"), float("nan"), float("nan"), -1)
    diff = ene_model[valid] - ene_ref[valid]
    imax = int(np.argmax(np.abs(diff)))
    return ErrorStatistics(
        n=n,
        rmse=float(np.sqrt(np.mean(diff**2))),
        mae=float(np.mean(np.abs(diff))),
        max_error=float(np.abs(diff[imax])),
        max_error_sample=int(samples[valid][imax]),
    )


def al_step_error_statistics(
    aldir: str, nreplicas: int
) -> Tuple[List[Tuple[int, ErrorStatistics]], ErrorStatistics]:
    """Error statistics of every replica of ``AL<n>`` and of all replicas together."""
    per_replica = []
    rows = []
    for replica in range(nreplicas):
        corr = read_energy_corr(os.path.join(aldir, str(replica), "energy_corr.dat"))
        per_replica.append((replica, error_statistics(corr)))
        rows.append(corr)
    overall = error_statistics(np.concatenate(rows) if rows else np.zeros((0, 3)))
    return per_replica, overall


def format_convergence(
    label: str,
    per_replica: Sequence[Tuple[int, ErrorStatistics]],
    overall: ErrorStatistics,
) -> List[str]:
    """Lines of ``convergence.dat`` (also suitable for logging)."""
    lines = [
        f"# {label}: errors of the model energies against the reference solver",
        "# (per supercell, in the units of energy_corr.dat)",
        "# replica  n  rmse  mae  max_error  max_error_sample",
    ]
    for replica, stats in per_replica:
        lines.append(_format_stats_row(str(replica), stats))
    lines.append(_format_stats_row("all", overall))
    return lines


def _format_stats_row(name: str, stats: ErrorStatistics) -> str:
    return (
        f"{name}  {stats.n}  {stats.rmse:.6g}  {stats.mae:.6g}  "
        f"{stats.max_error:.6g}  {stats.max_error_sample}"
    )


def write_convergence(path: str, lines: Iterable[str]) -> None:
    """Write the lines of :func:`format_convergence` to ``convergence.dat``."""
    with open(path, "w") as fi:
        for line in lines:
            fi.write(line + "\n")
