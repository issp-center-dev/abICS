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
Integrated execution of the reference-data step (``abics_mlref`` with
``[mlref] exec_mode = "integrated"``).

In the staged mode ``abics_mlref`` only writes solver inputs and parses
solver outputs; the solver itself is run by the user between two
invocations.  In the integrated mode a single, non-MPI ``abics_mlref``
process does everything for one active-learning step: it generates the
configurations, runs the reference solver for each of them through
libEnsemble's ``MPIExecutor`` (``nworkers`` runs at a time, each with
``nprocs`` MPI ranks), chains the relaxation stages, retries failed
runs, and writes the same files as the staged mode
(``AL<n>/<replica>/structure.<i>.vasp``, ``energy_corr.dat``,
``relax_max.dat``, ``convergence.dat``, ``ALloop.progress``).

The directory layout of the staged mode is kept: the solver runs in
``AL<n>/<replica>/input<i>/baseinput<k>`` (or in a copy under
``scratch_dir``, with the outputs copied back), so the results can be
used, inspected and resumed exactly like those of the staged mode.
A stage is complete when the marker file ``abics_mlref.done`` has been
written next to its output, which happens only after the solver exited
successfully and its output was parsed; re-running ``abics_mlref`` after
an interruption skips the complete stages and reruns the others.

libEnsemble is an optional dependency (``pip install libensemble``).
"""

from __future__ import annotations

import datetime
import fnmatch
import logging
import os
import shutil
import sys
import time
from dataclasses import asdict, dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import numpy.random as rand
from pymatgen.core import Structure

from abics.applications.latgas_abinitio_interface import activelearn as al
from abics.applications.latgas_abinitio_interface.base_solver import create_solver
from abics.applications.latgas_abinitio_interface.defect import DFTConfigParams, defect_config
from abics.applications.latgas_abinitio_interface.params import ALParams, DFTParams, TrainerParams
from abics.mlref import RefParams

logger = logging.getLogger("main")

SOLVER_APP = "reference_solver"
DONE_MARKER = "abics_mlref.done"
"""Written into ``baseinput<k>`` when the solver exited successfully and its output was parsed."""


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------


@dataclass
class IntegratedParams:
    """``[mlref.integrated]`` section."""

    nworkers: int = 1
    "Number of solver runs executed at the same time (libEnsemble workers)"
    nprocs: int = 1
    "MPI ranks per solver run"
    nthreads: int = 1
    "Threads per MPI rank (passed to the solver's command line generator)"
    timeout: float = 0.0
    "Seconds after which a solver run is killed (0: no limit)"
    retries: int = 1
    "How many times a failed (or timed-out) solver run is repeated"
    output_wait: float = 10.0
    "Seconds to keep re-trying to read the output after the solver exited (shared file systems)"
    scratch_dir: str = ""
    "Directory in which the solver runs (must be visible from the solver's nodes); empty: run in place"
    copy_back_exclude: List[str] = field(default_factory=list)
    "Glob patterns of files not copied back from scratch_dir (e.g. '*.wfc*')"
    cores: int = 0
    "Cores (libEnsemble resource sets) available to the solver runs (0: nworkers * nprocs)"

    @classmethod
    def from_dict(cls, d: dict) -> "IntegratedParams":
        params = cls()
        params.nworkers = int(d.get("nworkers", 1))
        params.nprocs = int(d.get("nprocs", 1))
        params.nthreads = int(d.get("nthreads", 1))
        params.timeout = float(d.get("timeout", 0.0))
        params.retries = int(d.get("retries", 1))
        params.output_wait = float(d.get("output_wait", 10.0))
        params.scratch_dir = os.path.expandvars(os.path.expanduser(str(d.get("scratch_dir", ""))))
        params.copy_back_exclude = list(d.get("copy_back_exclude", []))
        params.cores = int(d.get("cores", 0))
        if params.nworkers < 1 or params.nprocs < 1 or params.nthreads < 1:
            raise ValueError("[mlref.integrated]: nworkers, nprocs and nthreads must be >= 1")
        if params.retries < 0 or params.timeout < 0 or params.output_wait < 0 or params.cores < 0:
            raise ValueError("[mlref.integrated]: retries, timeout, output_wait and cores must be >= 0")
        if 0 < params.cores < params.nprocs:
            raise ValueError(
                f"[mlref.integrated]: cores = {params.cores} is smaller than nprocs = {params.nprocs}; "
                "no solver run could ever be scheduled"
            )
        if any(c.isspace() for c in params.scratch_dir):
            raise ValueError("[mlref.integrated]: scratch_dir must not contain whitespace "
                             "(libEnsemble splits the solver command line on whitespace)")
        return params


def require_libensemble() -> None:
    """Raise a clear error when libEnsemble is not installed."""
    try:
        import libensemble  # noqa: F401
    except ImportError as e:
        raise RuntimeError(
            "[mlref] exec_mode = 'integrated' requires libEnsemble "
            "(pip install libensemble); use exec_mode = 'staged' otherwise."
        ) from e


def resolve_solver_path(path: str) -> str:
    """Absolute path of the solver executable (``path`` may be a command on PATH)."""
    if os.path.isfile(path):
        return os.path.abspath(path)
    found = shutil.which(path)
    if found is None:
        raise FileNotFoundError(
            f"[mlref.solver] path = '{path}' is not an existing file nor a command on PATH"
        )
    return os.path.abspath(found)


MPI_LAUNCHER_VARIABLES = ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMI_RANK", "MPI_LOCALNRANKS", "MV2_COMM_WORLD_SIZE")


def launched_by_mpi() -> bool:
    """True when this process was started by mpiexec/mpirun/srun (even with one rank).

    Checked through the environment so that MPI is not initialized here.
    Running under a launcher is refused regardless of the rank count: the
    nested ``mpirun`` of the solver runs would fail or misbehave.
    """
    return any(var in os.environ for var in MPI_LAUNCHER_VARIABLES)


# ---------------------------------------------------------------------------
# One configuration: run the relaxation stages
# ---------------------------------------------------------------------------


@dataclass
class ReferenceTask:
    replica: int
    index: int
    "Sample index (the number in input<i>)"
    input_dir: str


@dataclass
class TaskResult:
    replica: int
    index: int
    status: str
    "'ok', or a description of the failure"
    energy: float = float("nan")
    dmax_id: int = -1
    dmax: float = float("nan")
    attempts: int = 0


Launcher = Callable[[str, int], str]
"""``launch(run_dir, attempt) -> status`` runs the solver in ``run_dir``;
``status`` is ``"ok"`` or a description of the failure."""


def try_parse(solver_output, stage_dir: str) -> Optional[Tuple[float, Structure]]:
    """Energy and structure of a finished solver run, or ``None`` if it cannot be parsed.

    Some solver parsers call ``exit()`` on a missing output, so ``SystemExit``
    is caught as well.
    """
    try:
        energy, st = solver_output.get_results(stage_dir)
    except (Exception, SystemExit):  # any parse error means "not finished"
        return None
    if energy is None or not np.isfinite(energy):
        return None
    return float(energy), st


def parse_with_wait(solver_output, stage_dir: str, output_wait: float) -> Optional[Tuple[float, Structure]]:
    """Parse the output, re-trying for ``output_wait`` seconds.

    On a shared file system the files written by the solver's node may
    become visible to this process a little later than the exit status.
    """
    deadline = time.time() + output_wait
    while True:
        parsed = try_parse(solver_output, stage_dir)
        if parsed is not None or time.time() >= deadline:
            return parsed
        time.sleep(min(1.0, max(0.1, deadline - time.time())))


def stage_done(stage_dir: str) -> bool:
    return os.path.isfile(os.path.join(stage_dir, DONE_MARKER))


def mark_stage_done(stage_dir: str, attempt: int) -> None:
    with open(os.path.join(stage_dir, DONE_MARKER), "w") as f:
        f.write(f"finished {datetime.datetime.now().isoformat()} after {attempt} attempt(s)\n")


def run_task_stages(
    task: ReferenceTask,
    solver,
    alparams: ALParams,
    retries: int,
    launch: Launcher,
    output_wait: float = 0.0,
) -> Tuple[str, float, Optional[Structure], int]:
    """Run all relaxation stages of one configuration.

    Returns ``(status, energy, structure, attempts)``; ``structure`` is the
    solver's final structure without the ignored species.  A stage whose
    marker file (:data:`DONE_MARKER`) exists is not run again (resume after
    interruption); a stage with output but no marker is rerun, because the
    output may be that of an interrupted run.
    """
    solver_input = solver.input
    solver_output = solver.output
    energy = float("nan")
    structure: Optional[Structure] = None
    attempts = 0
    nstages = len(alparams.base_input_dir)
    for stage in range(nstages):
        sdir = al.stage_dir(task.input_dir, stage)
        if stage > 0 and not os.path.isdir(sdir):
            solver_input.from_directory(alparams.base_input_dir[stage])
            al.write_next_stage_input(solver_input, structure, task.input_dir, stage)
        parsed = try_parse(solver_output, sdir) if stage_done(sdir) else None
        if parsed is None and stage_done(sdir):
            os.remove(os.path.join(sdir, DONE_MARKER))  # marker without readable output: rerun
        status = "ok"
        attempt = 0
        while parsed is None and attempt <= retries:
            attempt += 1
            attempts += 1
            status = launch(sdir, attempt)
            if status == "ok":
                parsed = parse_with_wait(solver_output, sdir, output_wait)
                if parsed is None:
                    status = "output could not be parsed"
                else:
                    mark_stage_done(sdir, attempt)
        if parsed is None:
            return f"stage {stage}: {status}", energy, structure, attempts
        energy, structure = parsed
    return "ok", energy, structure, attempts


def finish_task(
    task: ReferenceTask,
    structure: Structure,
    alparams: ALParams,
    al_step: int,
    replica_dir: str,
) -> Tuple[int, float]:
    """Write ``structure.<i>.vasp`` as the staged mode does and measure the relaxation."""
    st_in = Structure.from_file(os.path.join(task.input_dir, "initial.vasp"))
    if alparams.ignore_species:
        al.restore_ignored_species(structure, st_in, alparams.ignore_species)
    if al_step == 0:
        # the staged mode sorts the structures of AL0 (and only those)
        structure.sort(key=lambda site: site.species_string)
    structure.to(fmt="POSCAR", filename=os.path.join(replica_dir, f"structure.{task.index}.vasp"))
    return al.max_displacement(st_in, structure)


ALWAYS_EXCLUDE = ("machinefile_autogen_*",)
"""Files of the MPI launcher that are never copied back from the run directory."""


def _excluded(relpath: str, patterns: Sequence[str]) -> bool:
    """True if the relative path, its basename or any of its directory components matches a pattern."""
    parts = relpath.split(os.sep)
    return any(
        fnmatch.fnmatch(relpath, pat) or any(fnmatch.fnmatch(part, pat) for part in parts)
        for pat in patterns
    )


def copy_tree_filtered(src: str, dst: str, exclude: Sequence[str]) -> None:
    """Copy the files of ``src`` into ``dst`` except those matching ``exclude``.

    A pattern is matched against each file's path relative to ``src``, its
    name, and every directory name on its path, so ``"*.save"`` excludes the
    whole ``pwscf.save/`` directory.
    """
    patterns = list(exclude) + list(ALWAYS_EXCLUDE)
    for root, dirs, files in os.walk(src):
        rel = os.path.relpath(root, src)
        # do not descend into excluded directories
        dirs[:] = [
            d for d in dirs if not _excluded(d if rel == "." else os.path.join(rel, d), patterns)
        ]
        target = dst if rel == "." else os.path.join(dst, rel)
        os.makedirs(target, exist_ok=True)
        for name in files:
            relpath = name if rel == "." else os.path.join(rel, name)
            if _excluded(relpath, patterns):
                continue
            shutil.copy2(os.path.join(root, name), os.path.join(target, name))


def remove_launcher_files(run_dir: str) -> None:
    """Remove the MPI launcher's machinefiles that libEnsemble writes into the run directory."""
    for name in os.listdir(run_dir):
        if any(fnmatch.fnmatch(name, pat) for pat in ALWAYS_EXCLUDE):
            try:
                os.remove(os.path.join(run_dir, name))
            except OSError:
                pass


# ---------------------------------------------------------------------------
# libEnsemble generator and simulator functions
# ---------------------------------------------------------------------------


def task_list_gen(H, persis_info, gen_specs, libE_info):
    """Emit every configuration once; the whole list is known up front."""
    ntasks = gen_specs["user"]["ntasks"]
    H_o = np.zeros(ntasks, dtype=gen_specs["out"])
    H_o["task_index"] = np.arange(ntasks)
    H_o["resource_sets"] = gen_specs["user"]["nprocs"]
    return H_o, persis_info


def run_reference_task(H, persis_info, sim_specs, libE_info):
    """Run the solver for one configuration on a libEnsemble worker.

    Never raises: any problem is reported in the ``status`` field so that
    the manager can list the failed configurations.
    """
    from libensemble.executors.executor import Executor, TimeoutExpired

    user = sim_specs["user"]
    rootdir = user["rootdir"]
    task_index = int(H["task_index"][0])
    H_o = np.zeros(1, dtype=sim_specs["out"])
    H_o["task_index"] = task_index
    H_o["energy"] = np.nan
    H_o["dmax_id"] = -1
    H_o["dmax"] = np.nan
    try:
        os.chdir(rootdir)
        if rootdir not in sys.path:
            sys.path.insert(0, rootdir)
        # libEnsemble passes ``user`` through pydantic, which turns dataclasses
        # into dicts: rebuild the objects here
        task = ReferenceTask(**user["tasks"][task_index])
        ip = IntegratedParams(**user["integrated"])
        alparams = ALParams.from_dict(user["mlref_solver"])  # paths relative to rootdir
        replica_dir = user["replica_dirs"][task.replica]
        solver = create_solver(alparams.solver, alparams)
        exctr = Executor.executor

        def launch(stage_dir: str, attempt: int) -> str:
            run_dir = stage_dir
            if ip.scratch_dir:
                # under the per-run scratch root made by the manager (removed when a
                # rerun starts, so copies of an interrupted run do not pile up)
                run_dir = os.path.join(
                    user["scratch_root"],
                    f"{task.replica}_{task.index}_{os.path.basename(stage_dir)}_{attempt}",
                )
                shutil.rmtree(run_dir, ignore_errors=True)
                shutil.copytree(stage_dir, run_dir)
            args = solver.input.cl_args(ip.nprocs, ip.nthreads, run_dir)
            status = "ok"
            try:
                os.chdir(run_dir)
                libe_task = exctr.submit(
                    app_name=SOLVER_APP,
                    num_procs=ip.nprocs,
                    app_args=" ".join(str(a) for a in args),
                    stdout=os.path.join(run_dir, "stdout"),
                    stderr=os.path.join(run_dir, "stderr"),
                    wait_on_start=True,
                )
                try:
                    libe_task.wait(timeout=ip.timeout if ip.timeout > 0 else None)
                    if not libe_task.success:
                        status = f"solver exited with status {libe_task.errcode}"
                except TimeoutExpired:
                    libe_task.kill()
                    status = f"timeout after {ip.timeout} s"
            finally:
                os.chdir(rootdir)
                if ip.scratch_dir:
                    # copy back whatever is there (also after a failure, for inspection)
                    copy_tree_filtered(run_dir, stage_dir, ip.copy_back_exclude)
                    shutil.rmtree(run_dir, ignore_errors=True)
                else:
                    remove_launcher_files(run_dir)
            return status

        status, energy, structure, attempts = run_task_stages(
            task, solver, alparams, ip.retries, launch, ip.output_wait
        )
        H_o["attempts"] = attempts
        if status == "ok" and structure is not None:
            dmax_id, dmax = finish_task(task, structure, alparams, user["al_step"], replica_dir)
            H_o["energy"] = energy
            H_o["dmax_id"] = dmax_id
            H_o["dmax"] = dmax
    except (Exception, SystemExit) as e:  # noqa: BLE001 - reported, not raised
        detail = str(e)
        if isinstance(e, OSError) and e.filename:
            # keep the message short and the file recognizable
            try:
                detail = f"{e.strerror}: {os.path.relpath(e.filename, rootdir)}"
            except ValueError:
                detail = f"{e.strerror}: {e.filename}"
        status = f"internal error: {type(e).__name__}: {detail}"
        try:
            os.chdir(rootdir)
        except OSError:
            pass
    H_o["status"] = status[:255]
    return H_o, persis_info


def run_reference_calculations(
    tasks: Sequence[ReferenceTask],
    mlref_solver: dict,
    ip: IntegratedParams,
    rootdir: str,
    replica_dirs: Sequence[str],
    al_step: int,
    libe_dir: str,
    scratch_root: str = "",
) -> List[TaskResult]:
    """Run the reference solver for all ``tasks`` through libEnsemble and collect the results.

    ``mlref_solver`` is the ``[mlref.solver]`` section (the workers build
    their own ``ALParams`` and solver from it); ``replica_dirs[r]`` is
    ``AL<n>/<r>``; ``scratch_root`` is the directory under ``scratch_dir``
    for this run (ignored when ``scratch_dir`` is empty).
    """
    require_libensemble()
    from libensemble import Ensemble
    from libensemble.executors.mpi_executor import MPIExecutor
    from libensemble.specs import ExitCriteria, GenSpecs, LibeSpecs, SimSpecs

    alparams = ALParams.from_dict(mlref_solver)
    exctr = MPIExecutor()
    exctr.register_app(full_path=resolve_solver_path(alparams.path), app_name=SOLVER_APP)

    cores = ip.cores if ip.cores > 0 else ip.nworkers * ip.nprocs
    sim_specs = SimSpecs(
        sim_f=run_reference_task,
        inputs=["task_index"],
        outputs=[
            ("task_index", int),
            ("status", "U256"),
            ("energy", float),
            ("dmax_id", int),
            ("dmax", float),
            ("attempts", int),
        ],
        user={
            # plain dicts/lists only: pydantic converts dataclasses anyway
            "tasks": [asdict(t) for t in tasks],
            "mlref_solver": dict(mlref_solver),
            "integrated": asdict(ip),
            "rootdir": rootdir,
            "replica_dirs": list(replica_dirs),
            "al_step": al_step,
            "scratch_root": scratch_root,
        },
    )
    gen_specs = GenSpecs(
        gen_f=task_list_gen,
        outputs=[("task_index", int), ("resource_sets", int)],
        user={"ntasks": len(tasks), "nprocs": ip.nprocs},
    )
    libE_specs = LibeSpecs(
        comms="local",
        nworkers=ip.nworkers,
        gen_on_manager=True,
        # one resource set per core; each solver run takes nprocs of them
        num_resource_sets=cores,
        resource_info={"cores_on_node": (cores, cores)},
        sim_dirs_make=False,
        use_workflow_dir=True,
        workflow_dir_path=libe_dir,
        reuse_output_dir=True,
    )
    ensemble = Ensemble(
        sim_specs=sim_specs,
        gen_specs=gen_specs,
        libE_specs=libE_specs,
        exit_criteria=ExitCriteria(sim_max=len(tasks)),
        executor=exctr,
    )
    H, _, _ = ensemble.run()
    if not ensemble.is_manager:  # pragma: no cover - comms="local" always returns on the manager
        sys.exit(0)

    results = []
    for row in H:
        if not row["sim_ended"]:
            continue
        task = tasks[int(row["task_index"])]
        results.append(
            TaskResult(
                replica=task.replica,
                index=task.index,
                status=str(row["status"]),
                energy=float(row["energy"]),
                dmax_id=int(row["dmax_id"]),
                dmax=float(row["dmax"]),
                attempts=int(row["attempts"]),
            )
        )
    return results


# ---------------------------------------------------------------------------
# The integrated abics_mlref
# ---------------------------------------------------------------------------


def newest_mc_index(rootdir: str) -> int:
    index = 0
    while os.path.exists(os.path.join(rootdir, f"MC{index}")):
        index += 1
    return index


def seed_replica(seed: int, replica: int) -> None:
    """Seed the random number generator as RX_MPI_init does for rank ``replica``.

    With ``seed = 0`` every replica gets a fresh seed from the operating
    system, as every independently started rank of the staged mode does.
    """
    if seed > 0:
        rand.seed(seed + replica * 137)
    else:
        rand.seed()


def check_resumable(
    replica_dirs: Sequence[str],
    fmtstr_by_replica: Sequence[str],
    ndata: int,
    energy_ref: Optional[Sequence[Dict[int, float]]],
) -> Tuple[List[ReferenceTask], List[str]]:
    """Tasks of an existing ``AL<n>`` and the reasons why it cannot be resumed (empty if it can).

    Every replica directory must exist and hold exactly ``ndata`` input
    directories, and for ``AL<n>`` with ``n > 0`` every input index must be a
    sample of the preceding MC run (``energy_ref[r]``).
    """
    tasks: List[ReferenceTask] = []
    problems: List[str] = []
    for r, replica_dir in enumerate(replica_dirs):
        if not os.path.isdir(replica_dir):
            problems.append(f"{replica_dir} does not exist")
            continue
        indices = al.list_input_indices(replica_dir)
        if len(indices) != ndata:
            problems.append(f"{replica_dir} has {len(indices)} input directories, expected ndata = {ndata}")
        if energy_ref is not None:
            missing = [i for i in indices if i not in energy_ref[r]]
            if missing:
                problems.append(
                    f"{replica_dir}: inputs {missing} have no sample in the preceding MC run "
                    "(MC output or [sampling] sampler changed?)"
                )
        for index in indices:
            input_dir = os.path.join(replica_dir, fmtstr_by_replica[r].format(index))
            if not os.path.isdir(input_dir):
                problems.append(
                    f"{replica_dir}: input directory for sample {index} is not named "
                    f"{os.path.basename(input_dir)} (number of MC samples changed?)"
                )
                continue
            tasks.append(ReferenceTask(r, index, input_dir))
    return tasks, problems


def main_integrated(params_root) -> None:
    require_libensemble()
    if launched_by_mpi():
        logger.error("abics_mlref with exec_mode = 'integrated' must not be run under mpiexec/mpirun/srun:")
        logger.error("it is a single process that launches the solver runs itself.")
        sys.exit(1)

    rxparams = RefParams.from_dict(params_root["mlref"])
    ip = IntegratedParams.from_dict(params_root["mlref"].get("integrated", {}))
    alparams = ALParams.from_dict(params_root["mlref"]["solver"])
    mcparams = DFTParams.from_dict(params_root["sampling"]["solver"])
    configparams = DFTConfigParams.from_dict(params_root["config"])
    trainerparams = TrainerParams.from_dict(params_root["train"])
    sampler_type = params_root["sampling"].get("sampler", "RXMC")
    if alparams.solver_run_scheme == "function" or not alparams.path:
        logger.error("[mlref.solver] must give 'path' (the solver command) and run_scheme = 'subprocess'")
        logger.error("for exec_mode = 'integrated'.")
        sys.exit(1)
    resolve_solver_path(alparams.path)
    rootdir = os.getcwd()
    if any(c.isspace() for c in rootdir):
        logger.error(f"The working directory '{rootdir}' contains whitespace; libEnsemble splits the")
        logger.error("solver command line on whitespace, so exec_mode = 'integrated' cannot be used here.")
        sys.exit(1)

    ndata = rxparams.ndata
    nreplicas = rxparams.nreplicas
    energy_ref_empty = trainerparams.energy_ref
    mult_now = np.prod(configparams.supercell)
    solver = create_solver(alparams.solver, alparams)
    solver_input = solver.input

    nextMC_index = newest_mc_index(rootdir)
    ALstep = nextMC_index
    aldir = os.path.join(rootdir, f"AL{ALstep}")
    replica_dirs = [os.path.join(aldir, str(r)) for r in range(nreplicas)]
    scratch_root = ""
    if ip.scratch_dir:
        # one directory per AL step: a rerun removes the copies left by an interrupted run
        scratch_root = os.path.join(ip.scratch_dir, f"abics_mlref_AL{ALstep}")
        shutil.rmtree(scratch_root, ignore_errors=True)
        os.makedirs(scratch_root, exist_ok=True)
    progress_path = os.path.join(rootdir, al.PROGRESS_FILE)
    if ALstep == 0:
        if os.path.exists(progress_path):
            logger.error("It seems you've already run the first active learning step. You should train now.")
            sys.exit(1)
    else:
        MCdirname = f"MC{ALstep - 1}"
        last_li = al.last_progress(rootdir)
        if last_li != MCdirname:
            logger.error("The last action:")
            logger.error(f"  expected: {MCdirname}")
            logger.error(f"  actual:   {last_li}")
            logger.error("You shouldn't run activelearn now.")
            logger.error("Either train (abics_train) or MC (abics_sampling) first.")
            sys.exit(1)

    logger.info(f"-Running AL{ALstep} in integrated mode: {ip.nworkers} solver run(s) at a time, "
                f"{ip.nprocs} MPI rank(s) each")

    # --- candidates of the preceding MC run (AL1 and later) -----------------
    energy_ref: Optional[List[Dict[int, float]]] = None
    candidates_by_replica: List[List[al.ReferenceCandidate]] = []
    if ALstep == 0:
        fmtstr_by_replica = [al.input_dir_format(ndata)] * nreplicas
    else:
        energy_ref = []
        fmtstr_by_replica = []
        for r in range(nreplicas):
            candidates = al.reference_candidates(os.path.join(rootdir, MCdirname, str(r)), sampler_type)
            candidates_by_replica.append(candidates)
            energy_ref.append({c.index: c.energy + energy_ref_empty * mult_now for c in candidates})
            fmtstr_by_replica.append(al.input_dir_format(len(candidates)))
        too_few = [(r, len(c)) for r, c in enumerate(candidates_by_replica) if len(c) < ndata]
        if too_few:
            for r, n in too_few:
                logger.error(f"  {MCdirname}/{r} has only {n} samples, fewer than [mlref] ndata = {ndata}")
            logger.error("-Reduce ndata or sample more configurations (nsteps / sample_frequency).")
            sys.exit(1)

    # --- configurations and solver inputs ----------------------------------
    tasks: List[ReferenceTask] = []
    if os.path.isdir(aldir):
        logger.info(f"-{aldir} exists: resuming, stages with a {DONE_MARKER} marker are not rerun")
        tasks, problems = check_resumable(replica_dirs, fmtstr_by_replica, ndata, energy_ref)
        if problems:
            for problem in problems:
                logger.error(f"  {problem}")
            logger.error(f"-AL{ALstep} is incomplete and cannot be resumed. Remove the directory to "
                         "generate the configurations again.")
            sys.exit(1)
    else:
        logger.info(f"-Preparing {alparams.solver} input in AL{ALstep}/*/input*/baseinput0")
        solver_input.from_directory(alparams.base_input_dir[0])
        rundirs = []
        for r in range(nreplicas):
            seed_replica(rxparams.seed, r)
            os.makedirs(replica_dirs[r])
            if ALstep == 0:
                # A fresh DFTConfigParams per replica, as every staged rank has: the
                # lattice of its base structure caches its inverse on first use, and a
                # cached inverse changes the last digit of the fractional coordinates
                # written for the following replicas.
                config = defect_config(DFTConfigParams.from_dict(params_root["config"]))
                structures = enumerate(al.generate_random_structures(config, ndata))
            else:
                selected = al.select_candidates(candidates_by_replica[r], rxparams)
                structures = zip(
                    [c.index for c in selected],
                    al.load_candidate_structures(selected, mcparams.ignore_species),
                )
            for index, structure in structures:
                input_dir = os.path.join(replica_dirs[r], fmtstr_by_replica[r].format(index))
                rundirs.append(al.write_reference_input(solver_input, structure, alparams.perturb, input_dir))
                tasks.append(ReferenceTask(r, index, input_dir))
        al.write_rundirs(rootdir, rundirs)

    # --- solver runs -------------------------------------------------------
    logger.info(f"-Running {alparams.solver} for {len(tasks)} configurations")
    t0 = time.time()
    results = run_reference_calculations(
        tasks,
        params_root["mlref"]["solver"],
        ip,
        rootdir,
        replica_dirs,
        ALstep,
        os.path.join(aldir, "libE"),
        scratch_root,
    )
    if scratch_root:
        shutil.rmtree(scratch_root, ignore_errors=True)
    logger.info(f"--{len(results)} solver runs finished in {time.time() - t0:.1f} s")
    failed = [res for res in results if res.status != "ok"]
    if failed or len(results) != len(tasks):
        for res in failed:
            logger.error(f"  AL{ALstep}/{res.replica}/{fmtstr_by_replica[res.replica].format(res.index)}: "
                         f"{res.status} (after {res.attempts} attempt(s))")
        if len(results) != len(tasks):
            logger.error(f"  {len(tasks) - len(results)} configuration(s) returned no result")
        logger.error(f"-{len(failed) + len(tasks) - len(results)} of {len(tasks)} solver runs failed. "
                     "Fix the cause and run abics_mlref again: finished stages are not rerun.")
        sys.exit(1)

    # --- per-replica result files, as in the staged mode ------------------
    by_replica: List[List[TaskResult]] = [[] for _ in range(nreplicas)]
    for res in results:
        by_replica[res.replica].append(res)
    replica_rows = []
    for r in range(nreplicas):
        rows = sorted(by_replica[r], key=lambda res: res.index)
        if ALstep == 0:
            corr = [("NaN", res.energy, res.index) for res in rows]
        else:
            corr = [(energy_ref[r][res.index], res.energy, res.index) for res in rows]
            al.write_relax_max(
                os.path.join(replica_dirs[r], "relax_max.dat"),
                [(res.dmax_id, res.dmax) for res in rows],
            )
        al.write_energy_corr(os.path.join(replica_dirs[r], "energy_corr.dat"), corr)
        replica_rows.append([(float("nan") if e == "NaN" else e, ref, i) for e, ref, i in corr])
    if ALstep > 0:
        per_replica, overall = al.al_step_error_statistics(replica_rows)
        lines = al.format_convergence(f"AL{ALstep}", per_replica, overall)
        al.write_convergence(os.path.join(aldir, "convergence.dat"), lines)
        logger.info(f"-Errors of the model energies against {alparams.solver} results (AL{ALstep}/convergence.dat):")
        for line in lines:
            logger.info(f"  {line}")
        if al.write_parity_plot(os.path.join(aldir, "parity.png"), f"AL{ALstep}", replica_rows):
            logger.info(f"-Parity plot written to AL{ALstep}/parity.png")
        else:
            logger.info("-matplotlib is not available; no parity plot written")

    logger.info("-Writing ALloop.progress")
    al.record_progress(rootdir, f"AL{ALstep}")
    logger.info("--Done. Please run abics_train next.")
    logger.info(f"Exiting normally on {datetime.datetime.now()}.\n")
