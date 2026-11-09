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
interruption (finished stages are not repeated; a Monte Carlo run that
was interrupted is moved aside and repeated) and ``max_iterations`` can be
raised to extend a finished loop.  Before each stage it checks the
remaining wall time and stops cleanly, optionally running a resubmission
command, when the stage would not finish in time.

After each reference-data step the errors of the model are collected in
``al_summary.dat``; when the target number of iterations is reached a
description of the final model is written to ``al_driver/model_metadata.json``.
With ``[al] evaluate_final_model = true`` one more reference-data step is
run once exactly ``max_iterations`` Monte Carlo runs exist, so that the
final model itself is compared with the reference solver
(``AL<max_iterations>/convergence.dat``); its data are used for training
only if the loop is extended later.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import logging
import os
import re
import shlex
import shutil
import signal
import subprocess
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

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
STAGES = ("mlref", "train", "sampling")
TERMINATE_GRACE = 30.0
"""Seconds between SIGTERM and SIGKILL when a stage exceeds stage_timeout."""


class DriverError(RuntimeError):
    """The loop cannot continue; the reasons were logged."""


class DriverInterrupted(BaseException):
    """SIGTERM (e.g. the scheduler's time limit) or SIGINT reached the driver."""


def install_signal_handlers() -> None:
    """Turn SIGTERM into :class:`DriverInterrupted` so that cleanup code runs.

    Without this the default disposition kills the process at once: the lock
    file would stay and a running stage would be orphaned.
    """

    def handler(signum, frame):
        raise DriverInterrupted(signal.Signals(signum).name)

    signal.signal(signal.SIGTERM, handler)


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
    sampling_command: str = "mpiexec -n {nprocs} abics_sampling {input}"
    "Command of the sampling stage; {nprocs} = nreplicas * number of ensemble members, {nreplicas} = [sampling] nreplicas"
    stage_timeout: float = 0.0
    "Seconds after which a training or sampling child process is terminated (0: no limit)"
    walltime: float = 0.0
    "Wall time available to this run of abics_al in seconds (0: unlimited)"
    walltime_margin: float = 600.0
    "Seconds kept in reserve: a stage is not started if it may not finish within walltime - margin"
    stage_estimates: Dict[str, float] = field(default_factory=dict)
    "Expected duration of each stage (mlref, train, sampling) in seconds, used until the stage has run once"
    progress_wait: float = 30.0
    "Seconds to wait for the progress written by a child process to become visible (shared file systems)"
    resubmit_command: str = ""
    "Command run when abics_al stops because of the wall time (e.g. 'sbatch job.sh')"
    generate_train_input: bool = True
    "Generate the aenet training input with abics_gen_aenet_input if train.base_input_dir is missing"
    evaluate_final_model: bool = False
    "Once exactly max_iterations Monte Carlo runs exist, run one more reference-data step (without training) to evaluate the final model"

    @classmethod
    def from_dict(cls, d: dict) -> "DriverParams":
        params = cls()
        params.max_iterations = int(d.get("max_iterations", 1))
        params.evaluate_final_model = bool(d.get("evaluate_final_model", False))
        params.train_command = str(d.get("train_command", params.train_command))
        params.sampling_command = str(d.get("sampling_command", params.sampling_command))
        params.stage_timeout = float(d.get("stage_timeout", 0.0))
        params.walltime = float(d.get("walltime", 0.0))
        params.walltime_margin = float(d.get("walltime_margin", 600.0))
        params.stage_estimates = {str(k): float(v) for k, v in dict(d.get("stage_estimates", {})).items()}
        params.progress_wait = float(d.get("progress_wait", 30.0))
        params.resubmit_command = str(d.get("resubmit_command", ""))
        params.generate_train_input = bool(d.get("generate_train_input", True))
        if params.max_iterations < 1:
            raise ValueError("[al] max_iterations must be >= 1")
        if min(params.stage_timeout, params.walltime, params.walltime_margin, params.progress_wait) < 0:
            raise ValueError("[al] stage_timeout, walltime, walltime_margin and progress_wait must be >= 0")
        unknown = sorted(set(params.stage_estimates) - set(STAGES))
        if unknown or any(v < 0 for v in params.stage_estimates.values()):
            raise ValueError(f"[al] stage_estimates: keys must be among {STAGES} with values >= 0 (got {unknown})")
        return params


def expand_command(template: str, **values) -> List[str]:
    """Split a stage command into arguments, substituting the known placeholders.

    Only ``{name}`` for the given ``values`` is replaced (other braces are
    left as they are, so launcher options such as ``--cpu-bind=map_cpu:{0,1}``
    work), and ``~`` and environment variables are expanded in every
    argument.  No shell is involved.
    """
    args = []
    for token in shlex.split(template):
        token = os.path.expanduser(os.path.expandvars(token))
        for name, value in values.items():
            token = token.replace("{" + name + "}", str(value))  # within one argument: spaces survive
        args.append(token)
    return args


# ---------------------------------------------------------------------------
# State of the loop
# ---------------------------------------------------------------------------


def read_progress(rootdir: str) -> List[str]:
    path = os.path.join(rootdir, al.PROGRESS_FILE)
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def wait_for_progress(rootdir: str, expected_last: str, timeout: float) -> bool:
    """Wait until ``ALloop.progress`` ends with ``expected_last`` (shared file systems lag)."""
    deadline = time.time() + timeout
    while True:
        progress = read_progress(rootdir)
        if progress and progress[-1] == expected_last:
            return True
        if time.time() >= deadline:
            return False
        time.sleep(min(1.0, max(0.1, deadline - time.time())))


def completed_mc_runs(progress: List[str]) -> int:
    """Number of Monte Carlo runs recorded as finished in ``ALloop.progress``."""
    return sum(1 for p in progress if p.startswith("MC"))


def trained_steps(progress: List[str]) -> int:
    """Number of AL steps whose data the current model was trained on (0: no model yet).

    The model in ``baseinput`` is the one made by the last ``train`` entry,
    which follows ``AL<k>`` and used AL0 .. AL<k>.
    """
    last_train = max((i for i, p in enumerate(progress) if p == "train"), default=None)
    if last_train is None:
        return 0
    al_before = [p for p in progress[:last_train] if p.startswith("AL")]
    return int(al_before[-1][2:]) + 1 if al_before else 0


def next_action(progress: List[str], newest_mc: int, max_iterations: int) -> Tuple[str, int]:
    """What to run next, from ``ALloop.progress`` and the existing ``MC<n>`` directories.

    Returns ``(action, n)`` with action ``"mlref"`` (run AL<n>), ``"train"``
    (train after AL<n>), ``"sampling"`` (run MC<n>), ``"sampling_interrupted"``
    (MC<n> exists but was never recorded: a sampling run that did not finish,
    to be moved aside and repeated) or ``"done"`` (``n`` finished Monte Carlo
    runs are recorded and ``n >= max_iterations``; a reference-data or
    training stage beyond the target is left as it is).
    Raises :class:`DriverError` when the progress file and the directories
    disagree in any other way.
    """
    if not progress:
        if newest_mc != 0:
            raise DriverError(f"MC0..MC{newest_mc - 1} exist but ALloop.progress does not")
        if max_iterations <= 0:
            return "done", 0
        return "mlref", 0
    last = progress[-1]
    completed = completed_mc_runs(progress)
    if last.startswith("AL"):
        n = int(last[2:])
        if newest_mc != n:
            raise DriverError(f"ALloop.progress ends with {last} but {newest_mc} MC run(s) exist")
        if completed >= max_iterations:
            return "done", completed
        return "train", n
    if last == "train":
        al_entries = [p for p in progress if p.startswith("AL")]
        if not al_entries:
            raise DriverError("ALloop.progress ends with 'train' but has no AL entry")
        n = int(al_entries[-1][2:])
        if newest_mc == n:
            if completed >= max_iterations:
                return "done", completed
            return "sampling", n
        if newest_mc == n + 1:
            # abics_sampling creates MC<n> at its start and records it only at its end
            if completed >= max_iterations:
                return "done", completed
            return "sampling_interrupted", n
        raise DriverError(f"ALloop.progress ends with 'train' after AL{n} but {newest_mc} MC run(s) exist")
    if last.startswith("MC"):
        n = int(last[2:])
        if newest_mc != n + 1:
            raise DriverError(f"ALloop.progress ends with {last} but {newest_mc} MC run(s) exist")
        if completed >= max_iterations:
            return "done", completed
        return "mlref", n + 1
    raise DriverError(f"Unknown entry '{last}' in ALloop.progress")


RECENT_ACTIVITY = 120.0
"""Seconds: an unrecorded MC<n> modified more recently than this may still be written by a live run."""


def last_modified(path: str) -> float:
    """Newest modification time of ``path``, its replica directories and their files (one level)."""
    newest = os.path.getmtime(path)
    for entry in os.scandir(path):
        newest = max(newest, entry.stat().st_mtime)
        if entry.is_dir():
            for sub in os.scandir(entry.path):
                newest = max(newest, sub.stat().st_mtime)
    return newest


def set_aside_interrupted_mc(rootdir: str, n: int, driver_dir: str) -> str:
    """Move an unfinished ``MC<n>`` into the driver directory and return its new path.

    Refuses when the directory was modified within :data:`RECENT_ACTIVITY`
    seconds: the sampling run may still be alive although ``stage.pid``
    (:class:`StageRecord`) did not say so, for example MPI ranks that
    survived their killed launcher, a driver killed before the record was
    written, or directories from an abICS version without the record.
    """
    source = os.path.join(rootdir, f"MC{n}")
    age = time.time() - last_modified(source)
    if age < RECENT_ACTIVITY:
        raise DriverError(
            f"MC{n} is not recorded in ALloop.progress but was modified {age:.0f} s ago: a sampling run may "
            f"still be active. Make sure no abics_sampling is running, then rerun abics_al "
            f"(after {RECENT_ACTIVITY:.0f} s without changes the directory is moved aside and the sampling repeated)."
        )
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    target = os.path.join(driver_dir, f"MC{n}.interrupted-{stamp}")
    shutil.move(source, target)
    return target


LOCK_FILE = "abics_al.lock"


class DriverLock:
    """A lock file so that two drivers (or a driver and its own orphaned rerun) do not share a directory."""

    def __init__(self, driver_dir: str):
        self.path = os.path.join(driver_dir, LOCK_FILE)

    def acquire(self) -> None:
        if os.path.exists(self.path):
            try:
                with open(self.path) as f:
                    host, pid = f.read().split()[:2]
                pid = int(pid)
            except (OSError, ValueError):
                host, pid = "?", -1
            if host == os.uname().nodename and pid > 0 and _process_alive(pid):
                # the pid may have been reused by an unrelated process since the lock was written
                command = _process_command(pid)
                if not command or _is_driver_command(command):
                    raise DriverError(
                        f"another abics_al (pid {pid} on {host}) is running in this directory; "
                        f"stop it first or remove {self.path} if it is stale"
                    )
                logger.warning(f"-The lock file {self.path} names pid {pid}, which is now '{command}', not abics_al: "
                               "treating the lock as stale")
            elif host != os.uname().nodename:
                raise DriverError(
                    f"{self.path} says abics_al (pid {pid}) runs on {host}, which cannot be checked from here; "
                    "if that process is not running any more (e.g. its job was killed), remove the file and rerun"
                )
            else:
                logger.warning(f"-Removing the stale lock file {self.path} (pid {pid} is not running)")
        with open(self.path, "w") as f:
            f.write(f"{os.uname().nodename} {os.getpid()} {datetime.datetime.now().isoformat(timespec='seconds')}\n")

    def release(self) -> None:
        try:
            os.remove(self.path)
        except OSError:
            pass


def _process_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _process_command(pid: int) -> str:
    """Command line of a live process, or an empty string if it cannot be read."""
    try:
        out = subprocess.run(["ps", "-o", "command=", "-p", str(pid)], capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return ""
    return out.stdout.strip()


def _is_driver_command(command: str) -> bool:
    """True if a command line is an ``abics_al`` (the console script or the module)."""
    return re.search(r"(^|/)abics_al(\s|$)|abics\.scripts\.aldriver", command) is not None


def _group_alive(pgid: int) -> bool:
    """True while any process of the process group exists (also when it belongs to another user)."""
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _group_commands(pgid: int) -> List[str]:
    """Command lines of the processes of a process group (empty if none or if they cannot be listed)."""
    try:
        pids = subprocess.run(["pgrep", "-g", str(pgid)], capture_output=True, text=True, timeout=10).stdout.split()
    except (OSError, subprocess.SubprocessError):
        return []
    return [command for command in (_process_command(int(pid)) for pid in pids) if command]


STAGE_FILE = "stage.pid"


class StageRecord:
    """``al_driver/stage.pid``: the process group of the training or sampling stage that is running.

    The stage command is started in a session of its own, so its process
    group holds the command and, for an MPI launcher, the launcher process;
    the ranks started by ``mpiexec``/``srun`` are usually in groups of their
    own (or on other nodes) and are stopped by the launcher.  The record is
    written when a stage starts and removed when it ends, so that a rerun
    after the driver was killed without a chance to stop the stage (SIGKILL,
    a crash of the driver) can tell that the stage may still be running
    instead of starting it a second time.  Ranks that survive their launcher
    are not covered; the recent-modification rule of
    :func:`set_aside_interrupted_mc` remains the guard for them.
    """

    def __init__(self, driver_dir: str):
        self.path = os.path.join(driver_dir, STAGE_FILE)

    def write(self, stage: str, pgid: int, command: str = "") -> None:
        """Record the group; ``command`` is the stage's program (to recognize the group later)."""
        with open(self.path, "w") as f:
            f.write(f"{os.uname().nodename} {pgid} {stage} "
                    f"{datetime.datetime.now().isoformat(timespec='seconds')} {os.path.basename(command)}\n")

    def clear(self) -> None:
        try:
            os.remove(self.path)
        except OSError:
            pass

    def check(self) -> None:
        """Raise :class:`DriverError` if the recorded stage may still be running; remove a stale record."""
        if not os.path.exists(self.path):
            return
        try:
            with open(self.path) as f:
                fields = f.read().split()
            host, pgid, stage, started = fields[:4]
            command = fields[4] if len(fields) > 4 else ""
            pgid = int(pgid)
        except (OSError, ValueError, IndexError):
            logger.warning(f"-Removing the unreadable {self.path}")
            self.clear()
            return
        if host != os.uname().nodename:
            raise DriverError(
                f"{self.path} says the {stage} stage (process group {pgid}, started {started}) was running on {host}, "
                "which cannot be checked from here; make sure it has ended (e.g. its job was killed), then remove "
                f"the file (and {LOCK_FILE} next to it, if present) and rerun"
            )
        if pgid != os.getpgrp() and _group_alive(pgid):
            # a group id can be reused (after a reboot, or when pids wrapped around): look at
            # what runs in it before refusing
            commands = _group_commands(pgid)
            ours = [c for c in commands if "abics_" in c or (command and command.lower() in c.lower())]
            if ours or not commands:
                detail = f"; it runs: {ours[0]}" if ours else ""
                raise DriverError(
                    f"the {stage} stage started {started} by an earlier abics_al is still running (process group "
                    f"{pgid} on {host}{detail}); wait for it to end, or stop it with 'kill -TERM -- -{pgid}', then "
                    f"rerun. If no such stage exists any more (e.g. after a reboot), remove {self.path} and rerun"
                )
            logger.warning(f"-Process group {pgid} of {self.path} now runs '{commands[0]}', not the {stage} stage: "
                           "treating the record as stale")
        else:
            logger.warning(f"-Removing the stale {self.path}: the {stage} stage (process group {pgid}) started "
                           f"{started} is not running any more")
        self.clear()


# ---------------------------------------------------------------------------
# Wall time
# ---------------------------------------------------------------------------


class WallClock:
    """Remaining wall time of this run and the durations of earlier stages.

    The expected duration of a stage is the longest it has taken in any
    earlier run (stored in ``al_driver/durations.json``), or the user's
    ``stage_estimates`` until it has run once.
    """

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

    def expected(self, stage: str) -> Optional[float]:
        if stage in self.durations:
            return self.durations[stage]
        return self.params.stage_estimates.get(stage)

    def can_start(self, stage: str) -> Tuple[bool, str]:
        """Whether ``stage`` is expected to finish within the wall time.

        Returns ``(ok, message)``; the message is a warning when the stage
        is started although its duration is unknown.
        """
        if self.params.walltime <= 0:
            return True, ""
        remaining = self.params.walltime - self.elapsed()
        need = self.expected(stage)
        if need is None:
            if remaining <= self.params.walltime_margin:
                return False, f"only {remaining:.0f} s of wall time left (margin {self.params.walltime_margin:.0f} s)"
            return True, (
                f"duration of the {stage} stage unknown (no earlier run, no [al] stage_estimates): "
                f"starting it with {remaining:.0f} s left"
            )
        if need + self.params.walltime_margin > remaining:
            return False, (
                f"{remaining:.0f} s of wall time left, the {stage} stage is expected to take {need:.0f} s "
                f"and {self.params.walltime_margin:.0f} s are kept in reserve"
            )
        return True, ""

    def record(self, stage: str, seconds: float) -> None:
        self.durations[stage] = max(seconds, self.durations.get(stage, 0.0))
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        with open(self.path, "w") as f:
            json.dump(self.durations, f, indent=1)


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def run_child(
    args: Sequence[str], logfile: str, timeout: float, record: Optional[StageRecord] = None, stage: str = ""
) -> int:
    """Run a stage command, appending its output to ``logfile``; returns the exit status.

    The command runs in a session of its own; ``record`` (``al_driver/stage.pid``)
    names its process group while the command runs (see :class:`StageRecord`).
    On ``timeout`` the group is sent SIGTERM (``mpirun``/``srun`` then stop
    their ranks) and SIGKILL after :data:`TERMINATE_GRACE` seconds; the
    status is then -1.
    """
    logger.info(f"--Command: {' '.join(shlex.quote(a) for a in args)}")
    logger.info(f"--Output appended to: {logfile}")
    with open(logfile, "a") as fo:
        fo.write(f"=== {datetime.datetime.now().isoformat(timespec='seconds')} {' '.join(args)}\n")
        fo.flush()
        try:
            proc = subprocess.Popen(args, stdout=fo, stderr=subprocess.STDOUT, start_new_session=True)
        except OSError as e:
            logger.error(f"--Cannot start the command: {e}")
            return 127
        try:
            if record is not None:
                record.write(stage, proc.pid, args[0])
            return proc.wait(timeout=timeout or None)
        except subprocess.TimeoutExpired:
            logger.error(f"--No completion after {timeout} s: terminating")
            _stop_child(proc)
            return -1
        except BaseException:
            # the driver itself is being stopped (SIGTERM, Ctrl-C, an error) or the record
            # could not be written: do not leave the stage running, it would keep writing
            # into the directory
            logger.error("--abics_al is being stopped: terminating the stage")
            _stop_child(proc)
            raise
        finally:
            # keep the record if the stage was not reaped (a second signal during _stop_child):
            # the next run's check() will see whether the group is still alive
            if record is not None and proc.poll() is not None:
                record.clear()


def _signal_group(proc: subprocess.Popen, sig: int) -> None:
    """Send ``sig`` to the process group of ``proc`` (its session, see :func:`run_child`)."""
    try:
        os.killpg(proc.pid, sig)
    except ProcessLookupError:
        pass


def _stop_child(proc: subprocess.Popen) -> None:
    """SIGTERM to the stage's process group, then SIGKILL after :data:`TERMINATE_GRACE` seconds.

    An MPI launcher stops its ranks on SIGTERM; ranks of a launcher that had
    to be killed may survive (they are not in the group, see :class:`StageRecord`).
    """
    _signal_group(proc, signal.SIGTERM)
    try:
        proc.wait(timeout=TERMINATE_GRACE)
    except subprocess.TimeoutExpired:
        _signal_group(proc, signal.SIGKILL)
        proc.wait()


def ensure_train_input(params_root: dict) -> None:
    """Generate the aenet training input when ``train.base_input_dir`` does not exist yet."""
    train = params_root.get("train", {})
    if train.get("type", "aenet") != "aenet":
        return
    bases = train.get("base_input_dir", "./aenet_train_input")
    if not isinstance(bases, list):
        bases = [bases]
    from abics.scripts.gen_aenet_input import main_impl as gen_aenet_input

    for base in bases:
        base = os.path.expanduser(os.path.expandvars(base))
        if os.path.isdir(base) and os.listdir(base):
            continue
        logger.warning(f"-Training input {base} does not exist: generating the default template with "
                       "abics_gen_aenet_input (check it, or set [al] generate_train_input = false)")
        gen_aenet_input(params_root, base, False)


def sampling_nprocs(params_root: dict) -> Tuple[int, int]:
    """``(nreplicas, nprocs)`` of abics_sampling: one rank per replica and ensemble member."""
    sampling = params_root["sampling"]
    nreplicas = int(sampling["nreplicas"])
    base = sampling.get("solver", {}).get("base_input_dir", ["./baseinput"])
    nensemble = len(base) if isinstance(base, list) else 1
    return nreplicas, nreplicas * nensemble


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
    nreplicas, nprocs = sampling_nprocs(params_root)
    args = expand_command(params.train_command, input=tomlfile, nreplicas=nreplicas, nprocs=nprocs)
    logfile = os.path.join(driver_dir, f"train_{n}.log")
    rc = run_child(args, logfile, params.stage_timeout, StageRecord(driver_dir), "train")
    if rc != 0:
        raise DriverError(f"AL{n}: training failed (exit status {rc}); see {logfile}")
    if not wait_for_progress(rootdir, "train", params.progress_wait):
        raise DriverError(f"AL{n}: training ended but ALloop.progress does not record it; see {logfile}")


def run_sampling(params: DriverParams, params_root: dict, tomlfile: str, n: int, driver_dir: str, rootdir: str) -> None:
    logger.info(f"=== MC{n}: sampling ===")
    nreplicas, nprocs = sampling_nprocs(params_root)
    args = expand_command(params.sampling_command, input=tomlfile, nreplicas=nreplicas, nprocs=nprocs)
    logfile = os.path.join(driver_dir, f"sampling_{n}.log")
    rc = run_child(args, logfile, params.stage_timeout, StageRecord(driver_dir), "sampling")
    if rc != 0:
        raise DriverError(f"MC{n}: sampling failed (exit status {rc}); see {logfile}")
    if not wait_for_progress(rootdir, f"MC{n}", params.progress_wait) or not os.path.isdir(
        os.path.join(rootdir, f"MC{n}")
    ):
        raise DriverError(f"MC{n}: sampling ended but ALloop.progress or MC{n} is missing; see {logfile}")


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


def count_reference_data(rootdir: str, iterations: int) -> int:
    """Number of reference calculations in AL0 .. AL<iterations-1> (rows of energy_corr.dat)."""
    total = 0
    for n in range(iterations):
        r = 0
        while os.path.isdir(os.path.join(rootdir, f"AL{n}", str(r))):
            path = os.path.join(rootdir, f"AL{n}", str(r), "energy_corr.dat")
            if os.path.exists(path):
                total += len(np.loadtxt(path, ndmin=2))
            r += 1
    return total


def evaluation_step(rootdir: str, trained_on: int) -> Optional[str]:
    """``"AL<k>"`` if the model trained on AL0 .. AL<k-1> has been compared with the reference solver.

    That comparison is the reference-data step after the model's Monte Carlo
    run, ``AL<k>`` with ``k = trained_on``; it exists when the loop went on
    or when ``evaluate_final_model`` ran it.
    """
    if trained_on > 0 and os.path.exists(os.path.join(rootdir, f"AL{trained_on}", "convergence.dat")):
        return f"AL{trained_on}"
    return None


def write_model_metadata(
    rootdir: str,
    params_root: dict,
    tomlfile: str,
    iterations: int,
    driver_dir: str,
    trained_on: Optional[int] = None,
) -> str:
    """Describe the current model: provenance, scope and accuracy (a small model card).

    ``iterations`` is the number of completed Monte Carlo runs and
    ``trained_on`` the number of AL steps whose data the model in the
    sampling solver's ``base_input_dir`` was trained on (``iterations`` when
    the loop ended after a Monte Carlo run; more when it was stopped after
    a later reference-data step and training).
    """
    if trained_on is None:
        trained_on = iterations
    config = params_root.get("config", {})
    sampling = params_root.get("sampling", {})
    train = params_root.get("train", {})
    mlref = params_root.get("mlref", {})
    species = sorted(
        {g.get("name") for d in config.get("defect_structure", []) for g in d.get("groups", []) if g.get("name")}
        | {b.get("type") for b in config.get("base_structure", []) if b.get("type")}
    )
    # AL<trained_on>, when it exists, is the evaluation of the current model (the step
    # run by evaluate_final_model, or the next iteration's reference data)
    evaluated_in = evaluation_step(rootdir, trained_on)
    accuracy = {}
    for n in range(1, trained_on + 1 if evaluated_in else trained_on):
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
    if trained_on == 0:
        accuracy_note = "no model has been trained yet"
    else:
        accuracy_note = (
            "AL<n> compares the model trained after AL<n-1> (the one used in MC<n-1>) with the "
            "reference solver on configurations sampled in MC<n-1>; the current model, trained on "
            f"AL0 .. AL{trained_on - 1}, "
        )
        if evaluated_in:
            accuracy_note += f"is evaluated in {evaluated_in}"
        else:
            accuracy_note += "has not been evaluated (set [al] evaluate_final_model = true to add that step)"
    metadata = {
        "abics_version": __version__,
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "input_file": os.path.abspath(tomlfile),
        "input_sha256": toml_hash,
        "accuracy_note": accuracy_note,
        "provenance": {
            "iterations": iterations,
            "trained_on_steps": trained_on,
            "evaluated_in": evaluated_in,
            # rows of energy_corr.dat of the steps the model was trained on; abics_train may
            # drop configurations it cannot map, so this is an upper bound of the training-set size
            "reference_configurations": count_reference_data(rootdir, trained_on),
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


def resubmit(params: DriverParams) -> None:
    """Run the resubmission command; a failure is an error (the loop would otherwise silently end)."""
    logger.info(f"-Running the resubmission command: {params.resubmit_command}")
    try:
        proc = subprocess.run(expand_command(params.resubmit_command), check=False)
    except OSError as e:
        raise DriverError(f"the resubmission command could not be started: {e}") from e
    if proc.returncode != 0:
        raise DriverError(f"the resubmission command exited with status {proc.returncode}")


def run_driver(
    params_root: dict,
    tomlfile: str,
    max_iterations: Optional[int] = None,
    forget_durations: bool = False,
    evaluate_final_model: Optional[bool] = None,
) -> str:
    """Run the active-learning loop; returns ``"done"`` or ``"walltime"``.

    ``max_iterations`` and ``evaluate_final_model`` override the ``[al]``
    section when given.  Raises :class:`DriverError` when a stage fails or
    the directory is in an inconsistent state.
    """
    require_libensemble()
    if launched_by_mpi():
        raise DriverError("abics_al must not be run under mpiexec/mpirun/srun: it launches the stages itself")
    for section in ("sampling", "mlref", "train", "config"):
        if section not in params_root:
            raise DriverError(f"the input file has no [{section}] section")
    if params_root["sampling"].get("reload", False):
        raise DriverError(
            "[sampling] reload = true is not supported by abics_al: an interrupted Monte Carlo run is "
            "moved aside and repeated, while reload would continue the previous MC<n-1>. Set reload = false."
        )
    params = DriverParams.from_dict(params_root.get("al", {}))
    if max_iterations is not None:
        if max_iterations < 1:
            raise DriverError("--max-iterations must be >= 1")
        params.max_iterations = max_iterations
    if evaluate_final_model is not None:
        params.evaluate_final_model = evaluate_final_model
    rootdir = os.getcwd()
    driver_dir = os.path.join(rootdir, DRIVER_DIR)
    os.makedirs(driver_dir, exist_ok=True)
    if forget_durations:
        try:
            os.remove(os.path.join(driver_dir, DURATIONS_FILE))
            logger.info("-Recorded stage durations forgotten")
        except OSError:
            pass
    clock = WallClock(params, driver_dir)
    tomlfile = os.path.abspath(tomlfile)
    install_signal_handlers()
    lock = DriverLock(driver_dir)
    lock.acquire()
    try:
        StageRecord(driver_dir).check()  # a stage left running by a killed driver
        state = _run_loop(params, params_root, tomlfile, rootdir, driver_dir, clock)
    except DriverInterrupted as e:
        logger.error(f"-abics_al interrupted by {e}; the running stage was stopped. Rerun to continue.")
        raise DriverError(f"interrupted by {e}") from None
    finally:
        lock.release()
    if state == "walltime" and params.resubmit_command:
        resubmit(params)  # after the lock is released: the next driver may start at once
    return state


def _run_loop(params, params_root, tomlfile, rootdir, driver_dir, clock) -> str:
    logger.info(f"-Target: {params.max_iterations} Monte Carlo run(s) (MC0 .. MC{params.max_iterations - 1})")
    stages_finished_here = 0
    while True:
        progress = read_progress(rootdir)
        action, n = next_action(progress, newest_mc_index(rootdir), params.max_iterations)
        if action == "done":
            if progress and not progress[-1].startswith("MC"):
                if params.evaluate_final_model and progress[-1] == f"AL{n}":
                    logger.info(f"-{n} Monte Carlo run(s) are recorded and the evaluation step AL{n} exists: "
                                "nothing to run.")
                else:
                    logger.info(f"-{n} Monte Carlo run(s) are recorded, which meets the target; the stage after "
                                f"them ({progress[-1]}) is kept as it is.")
                break
            if not params.evaluate_final_model or n != params.max_iterations:
                # a lowered target (n > max_iterations) stops the loop without spending
                # reference calculations on an evaluation
                break
            # exactly max_iterations runs and the last entry is MC<n-1>: AL<n> does not exist yet,
            # and running it (without training afterwards) evaluates the model used in MC<n-1>
            action = "evaluate"
        stage = {"sampling_interrupted": "sampling", "evaluate": "mlref"}.get(action, action)
        ok, message = clock.can_start(stage)
        if not ok:
            pending = (
                f"the evaluation step AL{n} (the target of {n} Monte Carlo runs is met; only this mlref stage is pending)"
                if action == "evaluate"
                else f"the {stage} stage of step {n}"
            )
            if stages_finished_here == 0:
                raise DriverError(
                    f"{pending} does not fit into the wall time even at the start of the "
                    f"run ({message}); raise [al] walltime, lower the stage estimate, or remove "
                    f"{os.path.join(DRIVER_DIR, DURATIONS_FILE)} if the recorded duration is outdated. "
                    "Not resubmitting, since this run made no progress."
                )
            logger.info(f"-Stopping before {pending}: {message}.")
            logger.info("-Run abics_al again to continue; finished stages are not repeated.")
            return "walltime"
        if message:
            logger.warning(f"-{message}")
        if action == "sampling_interrupted":
            moved = set_aside_interrupted_mc(rootdir, n, driver_dir)
            logger.warning(f"-MC{n} exists but was never recorded in ALloop.progress: the sampling run did not finish.")
            logger.warning(f"--Moved it to {moved}; running the sampling again.")
            action = "sampling"
        t_stage = time.time()
        if action in ("mlref", "evaluate"):
            if action == "evaluate":
                logger.info(f"-The target is met: evaluating the current model (used in MC{n - 1}) with one more "
                            f"reference-data step, AL{n}, without training on it")
            run_mlref(params_root, n)
            if n > 0:
                for line in write_summary(rootdir):
                    if not line.startswith("#"):
                        logger.info(f"  {SUMMARY_FILE}: {line}")
        elif action == "train":
            run_train(params, params_root, tomlfile, n, driver_dir, rootdir)
        else:
            run_sampling(params, params_root, tomlfile, n, driver_dir, rootdir)
        clock.record(stage, time.time() - t_stage)
        stages_finished_here += 1
        logger.info(f"--{stage} of step {n} finished in {time.time() - t_stage:.0f} s")
        if action == "evaluate":
            break

    progress = read_progress(rootdir)
    completed = completed_mc_runs(progress)
    trained_on = trained_steps(progress)
    write_summary(rootdir)
    path = write_model_metadata(rootdir, params_root, tomlfile, completed, driver_dir, trained_on)
    evaluated_in = evaluation_step(rootdir, trained_on)
    logger.info(f"-{completed} Monte Carlo run(s) completed (target {params.max_iterations}); the current model "
                f"was trained on {trained_on} AL step(s); model description written to {path}")
    if evaluated_in:
        logger.info(f"-The current model is evaluated in {evaluated_in}/convergence.dat (last row of {SUMMARY_FILE}).")
    elif trained_on > 0:
        logger.info("-The current model itself has not been evaluated: set [al] evaluate_final_model = true "
                    "(or pass --evaluate-final-model) to run one more reference-data step for it.")
    logger.info(f"-Errors per step are in {SUMMARY_FILE}. Raise [al] max_iterations and rerun abics_al to continue.")
    return "done"
