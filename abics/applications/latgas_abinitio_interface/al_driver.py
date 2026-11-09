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
The active-learning driver ``abics_al``.

One non-MPI process runs the whole active-learning loop::

    AL0 (abics_mlref, integrated) -> train -> MC0 -> AL1 -> train -> MC1 -> ...

until ``[al] max_iterations`` Monte Carlo runs exist.  The reference-data
step runs in this process through libEnsemble
(:func:`~abics.applications.latgas_abinitio_interface.activelearn_integrated.run_integrated_step`);
training and sampling are started as child processes with the commands
given in the ``[al]`` section, so that the MPI launcher of the machine
(``mpiexec``, ``srun``, ...) can be used.

The state of the loop is the state of the directory: ``ALloop.progress``
and the ``MC<n>`` directories written by the three stage programs.  The
driver reads them to decide what to do next, so it can be rerun after an
interruption (finished stages are not repeated) and ``max_iterations``
can be raised to extend a finished loop.  Before each stage it checks
the remaining wall time and stops cleanly, optionally running a
resubmission command, when the stage would not finish in time.

After each reference-data step the errors of the model are appended to
``al_summary.dat``; when the target number of iterations is reached a
description of the final model is written to ``al_driver/model_metadata.json``.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import logging
import os
import shlex
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from abics import __version__
from abics.applications.latgas_abinitio_interface import activelearn as al
from abics.applications.latgas_abinitio_interface.activelearn_integrated import (
    IntegratedStepError,
    launched_by_mpi,
    newest_mc_index,
    require_libensemble,
    run_integrated_step,
)

logger = logging.getLogger("main")

DRIVER_DIR = "al_driver"
SUMMARY_FILE = "al_summary.dat"
DURATIONS_FILE = "durations.json"
METADATA_FILE = "model_metadata.json"


class DriverError(RuntimeError):
    """The loop cannot continue; the reasons were logged."""


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------


@dataclass
class DriverParams:
    """``[al]`` section."""

    max_iterations: int = 1
    "Number of Monte Carlo runs to complete (MC0 ... MC<max_iterations-1>)"
    train_command: str = "abics_train {input}"
    "Command of the training stage; {input} is the input file"
    sampling_command: str = "mpiexec -n {nreplicas} abics_sampling {input}"
    "Command of the sampling stage; {nreplicas} is [sampling] nreplicas"
    stage_timeout: float = 0.0
    "Seconds after which a training or sampling child process is killed (0: no limit)"
    walltime: float = 0.0
    "Wall time available to this run of abics_al in seconds (0: unlimited)"
    walltime_margin: float = 600.0
    "Seconds kept in reserve: a stage is not started if it may not finish within walltime - margin"
    resubmit_command: str = ""
    "Command run when abics_al stops because of the wall time (e.g. 'sbatch job.sh')"
    generate_train_input: bool = True
    "Generate the aenet training input with abics_gen_aenet_input if train.base_input_dir is missing"

    @classmethod
    def from_dict(cls, d: dict) -> "DriverParams":
        params = cls()
        params.max_iterations = int(d.get("max_iterations", 1))
        params.train_command = str(d.get("train_command", params.train_command))
        params.sampling_command = str(d.get("sampling_command", params.sampling_command))
        params.stage_timeout = float(d.get("stage_timeout", 0.0))
        params.walltime = float(d.get("walltime", 0.0))
        params.walltime_margin = float(d.get("walltime_margin", 600.0))
        params.resubmit_command = str(d.get("resubmit_command", ""))
        params.generate_train_input = bool(d.get("generate_train_input", True))
        if params.max_iterations < 1:
            raise ValueError("[al] max_iterations must be >= 1")
        if params.stage_timeout < 0 or params.walltime < 0 or params.walltime_margin < 0:
            raise ValueError("[al] stage_timeout, walltime and walltime_margin must be >= 0")
        return params


# ---------------------------------------------------------------------------
# State of the loop
# ---------------------------------------------------------------------------


def read_progress(rootdir: str) -> List[str]:
    path = os.path.join(rootdir, al.PROGRESS_FILE)
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def next_action(progress: List[str], newest_mc: int, max_iterations: int) -> Tuple[str, int]:
    """What to run next, from ``ALloop.progress`` and the existing ``MC<n>`` directories.

    Returns ``(action, n)`` with action ``"mlref"`` (run AL<n>), ``"train"``
    (train after AL<n>), ``"sampling"`` (run MC<n>) or ``"done"`` (n Monte
    Carlo runs exist, nothing left to do).  Raises :class:`DriverError` when
    the progress file and the directories disagree.
    """
    if not progress:
        if newest_mc != 0:
            raise DriverError(f"MC0..MC{newest_mc - 1} exist but ALloop.progress does not")
        return "mlref", 0
    last = progress[-1]
    if last.startswith("AL"):
        n = int(last[2:])
        if newest_mc != n:
            raise DriverError(f"ALloop.progress ends with {last} but {newest_mc} MC run(s) exist")
        return "train", n
    if last == "train":
        n = newest_mc
        al_entries = [p for p in progress if p.startswith("AL")]
        if not al_entries or int(al_entries[-1][2:]) != n:
            raise DriverError(f"ALloop.progress ends with 'train' but AL{n} is not the last AL entry")
        return "sampling", n
    if last.startswith("MC"):
        n = int(last[2:])
        if newest_mc != n + 1:
            raise DriverError(f"ALloop.progress ends with {last} but {newest_mc} MC run(s) exist")
        if n + 1 >= max_iterations:
            return "done", n + 1
        return "mlref", n + 1
    raise DriverError(f"Unknown entry '{last}' in ALloop.progress")


# ---------------------------------------------------------------------------
# Wall time
# ---------------------------------------------------------------------------


class WallClock:
    """Remaining wall time of this run and the durations of earlier stages."""

    def __init__(self, params: DriverParams, driver_dir: str):
        self.params = params
        self.t0 = time.time()
        self.path = os.path.join(driver_dir, DURATIONS_FILE)
        self.durations: Dict[str, float] = {}
        if os.path.exists(self.path):
            try:
                with open(self.path) as f:
                    self.durations = {k: float(v) for k, v in json.load(f).items()}
            except (OSError, ValueError):
                self.durations = {}

    def elapsed(self) -> float:
        return time.time() - self.t0

    def can_start(self, stage: str) -> Tuple[bool, str]:
        """Whether ``stage`` is expected to finish within the wall time."""
        if self.params.walltime <= 0:
            return True, ""
        need = self.durations.get(stage, 0.0)
        remaining = self.params.walltime - self.elapsed()
        if need + self.params.walltime_margin > remaining:
            return False, (
                f"{remaining:.0f} s of wall time left, the last {stage} stage took {need:.0f} s "
                f"and {self.params.walltime_margin:.0f} s are kept in reserve"
            )
        return True, ""

    def record(self, stage: str, seconds: float) -> None:
        self.durations[stage] = seconds
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        with open(self.path, "w") as f:
            json.dump(self.durations, f, indent=1)


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def run_child(command: str, logfile: str, timeout: float) -> int:
    """Run a stage command with its output in ``logfile``; returns the exit status."""
    args = shlex.split(command)
    logger.info(f"--Command: {command}")
    logger.info(f"--Output: {logfile}")
    with open(logfile, "w") as fo:
        try:
            proc = subprocess.run(args, stdout=fo, stderr=subprocess.STDOUT, timeout=timeout or None)
        except subprocess.TimeoutExpired:
            logger.error(f"--Killed after {timeout} s")
            return -1
        except FileNotFoundError as e:
            logger.error(f"--Command not found: {e.filename}")
            return 127
    return proc.returncode


def ensure_train_input(params_root: dict) -> None:
    """Generate the aenet training input when ``train.base_input_dir`` does not exist yet."""
    train = params_root.get("train", {})
    if train.get("type", "aenet") != "aenet":
        return
    base = train.get("base_input_dir", "./aenet_train_input")
    if isinstance(base, list):
        base = base[0]
    if os.path.isdir(base) and os.listdir(base):
        return
    from abics.scripts.gen_aenet_input import main_impl as gen_aenet_input

    logger.info(f"-Training input {base} does not exist: generating it with abics_gen_aenet_input")
    gen_aenet_input(params_root, None, False)


def run_mlref(params_root: dict, n: int) -> None:
    logger.info(f"=== AL{n}: reference calculations (integrated abics_mlref) ===")
    params = dict(params_root)
    params["mlref"] = dict(params_root["mlref"])
    params["mlref"]["exec_mode"] = "integrated"
    try:
        step = run_integrated_step(params)
    except IntegratedStepError as e:
        raise DriverError(f"AL{n}: the reference-data step failed") from e
    if step != n:
        raise DriverError(f"abics_mlref ran AL{step}, AL{n} was expected")


def run_train(params: DriverParams, params_root: dict, tomlfile: str, n: int, driver_dir: str, rootdir: str) -> None:
    logger.info(f"=== AL{n}: training ===")
    if params.generate_train_input:
        ensure_train_input(params_root)
    command = params.train_command.format(input=tomlfile, nreplicas=params_root["sampling"]["nreplicas"])
    rc = run_child(command, os.path.join(driver_dir, f"train_{n}.log"), params.stage_timeout)
    progress = read_progress(rootdir)
    if rc != 0 or not progress or progress[-1] != "train":
        raise DriverError(f"AL{n}: training failed (exit status {rc}); see {driver_dir}/train_{n}.log")


def run_sampling(params: DriverParams, params_root: dict, tomlfile: str, n: int, driver_dir: str, rootdir: str) -> None:
    logger.info(f"=== MC{n}: sampling ===")
    nreplicas = params_root["sampling"]["nreplicas"]
    command = params.sampling_command.format(input=tomlfile, nreplicas=nreplicas)
    rc = run_child(command, os.path.join(driver_dir, f"sampling_{n}.log"), params.stage_timeout)
    progress = read_progress(rootdir)
    if rc != 0 or not progress or progress[-1] != f"MC{n}" or not os.path.isdir(os.path.join(rootdir, f"MC{n}")):
        raise DriverError(f"MC{n}: sampling failed (exit status {rc}); see {driver_dir}/sampling_{n}.log")


# ---------------------------------------------------------------------------
# Summary and model description
# ---------------------------------------------------------------------------


def read_convergence(path: str) -> Dict[str, List[str]]:
    """Rows of ``convergence.dat`` keyed by the replica column (``"all"`` for the total)."""
    rows = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            rows[fields[0]] = fields[1:]
    return rows


def write_summary(rootdir: str) -> List[str]:
    """Write ``al_summary.dat``: one row per AL step with the errors over all replicas."""
    lines = [
        "# errors of the model energies against the reference solver, all replicas",
        "# (from AL<n>/convergence.dat; AL0 has no model yet)",
        "# step  n  rmse  mae  max_error",
    ]
    n = 1
    while os.path.isdir(os.path.join(rootdir, f"AL{n}")):
        path = os.path.join(rootdir, f"AL{n}", "convergence.dat")
        if os.path.exists(path):
            total = read_convergence(path).get("all")
            if total:
                lines.append(f"{n}  {total[0]}  {total[1]}  {total[2]}  {total[3]}")
        n += 1
    with open(os.path.join(rootdir, SUMMARY_FILE), "w") as f:
        f.write("\n".join(lines) + "\n")
    return lines


def count_training_data(rootdir: str) -> int:
    total = 0
    n = 0
    while os.path.isdir(os.path.join(rootdir, f"AL{n}")):
        r = 0
        while os.path.isdir(os.path.join(rootdir, f"AL{n}", str(r))):
            path = os.path.join(rootdir, f"AL{n}", str(r), "energy_corr.dat")
            if os.path.exists(path):
                total += len(np.loadtxt(path, ndmin=2))
            r += 1
        n += 1
    return total


def write_model_metadata(rootdir: str, params_root: dict, tomlfile: str, iterations: int, driver_dir: str) -> str:
    """Describe the final model: provenance, scope and accuracy (a small model card)."""
    config = params_root.get("config", {})
    sampling = params_root.get("sampling", {})
    train = params_root.get("train", {})
    mlref = params_root.get("mlref", {})
    species = sorted(
        {g.get("name") for d in config.get("defect_structure", []) for g in d.get("groups", []) if g.get("name")}
        | {b.get("type") for b in config.get("base_structure", []) if b.get("type")}
    )
    accuracy = {}
    for n in range(1, iterations):
        path = os.path.join(rootdir, f"AL{n}", "convergence.dat")
        if os.path.exists(path):
            total = read_convergence(path).get("all")
            if total:
                accuracy[f"AL{n}"] = {
                    "n": int(total[0]),
                    "rmse": float(total[1]),
                    "mae": float(total[2]),
                    "max_error": float(total[3]),
                    "units": "eV per supercell",
                }
    with open(tomlfile, "rb") as f:
        toml_hash = hashlib.sha256(f.read()).hexdigest()
    metadata = {
        "abics_version": __version__,
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "input_file": os.path.abspath(tomlfile),
        "input_sha256": toml_hash,
        "provenance": {
            "iterations": iterations,
            "training_configurations": count_training_data(rootdir),
            "reference_solver": {
                "type": mlref.get("solver", {}).get("type"),
                "base_input_dir": mlref.get("solver", {}).get("base_input_dir"),
                "perturb": mlref.get("solver", {}).get("perturb"),
            },
            "mlref_seed": mlref.get("seed", 0),
            "trainer": {
                "type": train.get("type"),
                "base_input_dir": train.get("base_input_dir"),
                "ignore_species": train.get("ignore_species"),
            },
            # lineage of a model fine-tuned from a pretrained one, when the trainer supports it
            "base_model": train.get("base_model"),
        },
        "scope": {
            "species": species,
            "ignored_species": sampling.get("solver", {}).get("ignore_species"),
            "supercell": config.get("supercell"),
            "kT_range": [sampling.get("kTstart"), sampling.get("kTend")],
            "nreplicas": sampling.get("nreplicas"),
            "sampler": sampling.get("sampler", "RXMC"),
        },
        "model_files": sampling.get("solver", {}).get("base_input_dir"),
        "accuracy": accuracy,
    }
    path = os.path.join(driver_dir, METADATA_FILE)
    with open(path, "w") as f:
        json.dump(metadata, f, indent=1)
    return path


# ---------------------------------------------------------------------------
# The loop
# ---------------------------------------------------------------------------


def run_driver(params_root: dict, tomlfile: str, max_iterations: Optional[int] = None) -> str:
    """Run the active-learning loop; returns ``"done"`` or ``"walltime"``.

    Raises :class:`DriverError` when a stage fails or the directory is in an
    inconsistent state.
    """
    require_libensemble()
    if launched_by_mpi():
        raise DriverError("abics_al must not be run under mpiexec/mpirun/srun: it launches the stages itself")
    params = DriverParams.from_dict(params_root.get("al", {}))
    if max_iterations is not None:
        params.max_iterations = max_iterations
    rootdir = os.getcwd()
    driver_dir = os.path.join(rootdir, DRIVER_DIR)
    os.makedirs(driver_dir, exist_ok=True)
    clock = WallClock(params, driver_dir)
    tomlfile = os.path.abspath(tomlfile)

    logger.info(f"-Target: {params.max_iterations} Monte Carlo run(s) (MC0 .. MC{params.max_iterations - 1})")
    while True:
        progress = read_progress(rootdir)
        action, n = next_action(progress, newest_mc_index(rootdir), params.max_iterations)
        if action == "done":
            break
        ok, reason = clock.can_start(action)
        if not ok:
            logger.info(f"-Stopping before {action} of step {n}: {reason}.")
            logger.info("-Run abics_al again to continue; finished stages are not repeated.")
            if params.resubmit_command:
                logger.info(f"-Running the resubmission command: {params.resubmit_command}")
                subprocess.run(shlex.split(params.resubmit_command), check=False)
            return "walltime"
        t_stage = time.time()
        if action == "mlref":
            run_mlref(params_root, n)
            if n > 0:
                for line in write_summary(rootdir):
                    if not line.startswith("#"):
                        logger.info(f"  {SUMMARY_FILE}: {line}")
        elif action == "train":
            run_train(params, params_root, tomlfile, n, driver_dir, rootdir)
        else:
            run_sampling(params, params_root, tomlfile, n, driver_dir, rootdir)
        clock.record(action, time.time() - t_stage)
        logger.info(f"--{action} of step {n} finished in {time.time() - t_stage:.0f} s")

    write_summary(rootdir)
    path = write_model_metadata(rootdir, params_root, tomlfile, params.max_iterations, driver_dir)
    logger.info(f"-{params.max_iterations} Monte Carlo run(s) completed; model description written to {path}")
    logger.info(f"-Errors per step are in {SUMMARY_FILE}. Raise [al] max_iterations and rerun abics_al to continue.")
    return "done"
