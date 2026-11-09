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
Re-running ``abics_mlref`` after an interruption skips the
configurations whose final solver output already parses.

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
    scratch_dir: str = ""
    "Directory (e.g. node-local) in which the solver runs; empty: run in place"
    copy_back_exclude: List[str] = field(default_factory=list)
    "Glob patterns of files not copied back from scratch_dir (e.g. '*.wfc*')"
    cores: int = 0
    "Cores available to the solver runs (0: nworkers * nprocs)"

    @classmethod
    def from_dict(cls, d: dict) -> "IntegratedParams":
        params = cls()
        params.nworkers = int(d.get("nworkers", 1))
        params.nprocs = int(d.get("nprocs", 1))
        params.nthreads = int(d.get("nthreads", 1))
        params.timeout = float(d.get("timeout", 0.0))
        params.retries = int(d.get("retries", 1))
        params.scratch_dir = str(d.get("scratch_dir", ""))
        params.copy_back_exclude = list(d.get("copy_back_exclude", []))
        params.cores = int(d.get("cores", 0))
        if params.nworkers < 1 or params.nprocs < 1 or params.retries < 0:
            raise ValueError("[mlref.integrated]: nworkers and nprocs must be >= 1 and retries >= 0")
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
    """Energy and structure of a finished solver run, or ``None`` if it cannot be parsed."""
    try:
        energy, st = solver_output.get_results(stage_dir)
    except Exception:  # noqa: BLE001 - any parse error means "not finished"
        return None
    if energy is None or not np.isfinite(energy):
        return None
    return float(energy), st


def run_task_stages(
    task: ReferenceTask,
    solver,
    alparams: ALParams,
    retries: int,
    launch: Launcher,
) -> Tuple[str, float, Optional[Structure], int]:
    """Run all relaxation stages of one configuration.

    Returns ``(status, energy, structure, attempts)``; ``structure`` is the
    solver's final structure without the ignored species.  Stages whose
    output already parses are not run again (resume after interruption).
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
        parsed = try_parse(solver_output, sdir)
        status = "ok"
        attempt = 0
        while parsed is None and attempt <= retries:
            attempt += 1
            attempts += 1
            status = launch(sdir, attempt)
            if status == "ok":
                parsed = try_parse(solver_output, sdir)
                if parsed is None:
                    status = "output could not be parsed"
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


def copy_tree_filtered(src: str, dst: str, exclude: Sequence[str]) -> None:
    """Copy the files of ``src`` into ``dst`` except those matching ``exclude`` patterns."""
    for root, dirs, files in os.walk(src):
        rel = os.path.relpath(root, src)
        target = dst if rel == "." else os.path.join(dst, rel)
        os.makedirs(target, exist_ok=True)
        for name in files:
            relpath = name if rel == "." else os.path.join(rel, name)
            if any(fnmatch.fnmatch(relpath, pat) or fnmatch.fnmatch(name, pat) for pat in exclude):
                continue
            shutil.copy2(os.path.join(root, name), os.path.join(target, name))


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
    """Run the solver for one configuration on a libEnsemble worker."""
    from libensemble.executors.executor import Executor, TimeoutExpired

    user = sim_specs["user"]
    rootdir = user["rootdir"]
    os.chdir(rootdir)
    if rootdir not in sys.path:
        sys.path.insert(0, rootdir)
    # libEnsemble passes ``user`` through pydantic, which turns dataclasses
    # into dicts: rebuild the objects here
    task = ReferenceTask(**user["tasks"][int(H["task_index"][0])])
    ip = IntegratedParams(**user["integrated"])
    alparams = ALParams.from_dict(user["mlref_solver"])  # paths relative to rootdir
    replica_dir = user["replica_dirs"][task.replica]
    solver = create_solver(alparams.solver, alparams)
    exctr = Executor.executor

    def launch(stage_dir: str, attempt: int) -> str:
        run_dir = stage_dir
        if ip.scratch_dir:
            run_dir = os.path.join(
                ip.scratch_dir,
                f"abics_{os.getpid()}_{task.replica}_{task.index}_{os.path.basename(stage_dir)}_{attempt}",
            )
            shutil.rmtree(run_dir, ignore_errors=True)
            shutil.copytree(stage_dir, run_dir)
        args = solver.input.cl_args(ip.nprocs, ip.nthreads, run_dir)
        os.chdir(run_dir)
        try:
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
            except TimeoutExpired:
                libe_task.kill()
                return f"timeout after {ip.timeout} s"
        finally:
            os.chdir(rootdir)
        if ip.scratch_dir:
            copy_tree_filtered(run_dir, stage_dir, ip.copy_back_exclude)
            shutil.rmtree(run_dir, ignore_errors=True)
        if libe_task.success:
            return "ok"
        return f"solver exited with status {libe_task.errcode}"

    status, energy, structure, attempts = run_task_stages(task, solver, alparams, ip.retries, launch)
    H_o = np.zeros(1, dtype=sim_specs["out"])
    H_o["task_index"] = int(H["task_index"][0])
    H_o["attempts"] = attempts
    if status == "ok" and structure is not None:
        dmax_id, dmax = finish_task(task, structure, alparams, user["al_step"], replica_dir)
        H_o["energy"] = energy
        H_o["dmax_id"] = dmax_id
        H_o["dmax"] = dmax
    else:
        H_o["energy"] = np.nan
        H_o["dmax_id"] = -1
        H_o["dmax"] = np.nan
    H_o["status"] = status
    return H_o, persis_info


def run_reference_calculations(
    tasks: Sequence[ReferenceTask],
    mlref_solver: dict,
    ip: IntegratedParams,
    rootdir: str,
    replica_dirs: Sequence[str],
    al_step: int,
    libe_dir: str,
) -> List[TaskResult]:
    """Run the reference solver for all ``tasks`` through libEnsemble and collect the results.

    ``mlref_solver`` is the ``[mlref.solver]`` section (the workers build
    their own ``ALParams`` and solver from it); ``replica_dirs[r]`` is
    ``AL<n>/<r>``.
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
            ("status", "U128"),
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
    """Seed the random number generator as RX_MPI_init does for rank ``replica``."""
    if seed > 0:
        rand.seed(seed + replica * 137)
    else:
        rand.seed(rand.randint(10000))


def main_integrated(params_root) -> None:
    require_libensemble()
    from mpi4py import MPI

    if MPI.COMM_WORLD.Get_size() > 1:
        logger.error("abics_mlref with exec_mode = 'integrated' must not be run under mpiexec:")
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

    ndata = rxparams.ndata
    nreplicas = rxparams.nreplicas
    energy_ref_empty = trainerparams.energy_ref
    mult_now = np.prod(configparams.supercell)
    solver = create_solver(alparams.solver, alparams)
    solver_input = solver.input
    rootdir = os.getcwd()

    nextMC_index = newest_mc_index(rootdir)
    ALstep = nextMC_index
    aldir = os.path.join(rootdir, f"AL{ALstep}")
    replica_dirs = {r: os.path.join(aldir, str(r)) for r in range(nreplicas)}
    progress_path = os.path.join(rootdir, al.PROGRESS_FILE)
    if ALstep == 0:
        if os.path.exists(progress_path):
            logger.error("It seems you've already run the first active learning step. You should train now.")
            sys.exit(1)
        fmtstr = al.input_dir_format(ndata)
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

    # --- configurations and solver inputs -------------------------------
    tasks: List[ReferenceTask] = []
    energy_ref: Dict[int, Dict[int, float]] = {}
    resume = os.path.isdir(aldir)
    if resume:
        logger.info(f"-{aldir} exists: resuming, configurations with finished solver runs are skipped")
    if ALstep > 0:
        for r in range(nreplicas):
            candidates = al.reference_candidates(os.path.join(rootdir, MCdirname, str(r)), sampler_type)
            energy_ref[r] = {c.index: c.energy + energy_ref_empty * mult_now for c in candidates}
            if r == 0:
                fmtstr = al.input_dir_format(len(candidates))
    if not resume:
        logger.info(f"-Preparing {alparams.solver} input in AL{ALstep}/*/input*/baseinput0")
        solver_input.from_directory(alparams.base_input_dir[0])
        rundirs = []
        for r in range(nreplicas):
            seed_replica(rxparams.seed, r)
            os.makedirs(replica_dirs[r])
            if ALstep == 0:
                config = defect_config(configparams)
                structures = enumerate(al.generate_random_structures(config, ndata))
            else:
                candidates = al.reference_candidates(os.path.join(rootdir, MCdirname, str(r)), sampler_type)
                selected = al.select_candidates(candidates, rxparams)
                structures = zip(
                    [c.index for c in selected],
                    al.load_candidate_structures(selected, mcparams.ignore_species),
                )
            for index, structure in structures:
                input_dir = os.path.join(replica_dirs[r], fmtstr.format(index))
                rundirs.append(al.write_reference_input(solver_input, structure, alparams.perturb, input_dir))
                tasks.append(ReferenceTask(r, index, input_dir))
        al.write_rundirs(rootdir, rundirs)
    else:
        for r in range(nreplicas):
            for index in al.list_input_indices(replica_dirs[r]):
                tasks.append(ReferenceTask(r, index, os.path.join(replica_dirs[r], fmtstr.format(index))))

    # --- solver runs -------------------------------------------------------
    logger.info(f"-Running {alparams.solver} for {len(tasks)} configurations")
    t0 = time.time()
    results = run_reference_calculations(
        tasks,
        params_root["mlref"]["solver"],
        ip,
        rootdir,
        [replica_dirs[r] for r in range(nreplicas)],
        ALstep,
        os.path.join(aldir, "libE"),
    )
    logger.info(f"--{len(results)} solver runs finished in {time.time() - t0:.1f} s")
    failed = [res for res in results if res.status != "ok"]
    if failed or len(results) != len(tasks):
        for res in failed:
            logger.error(f"  AL{ALstep}/{res.replica}/{fmtstr.format(res.index)}: {res.status} "
                         f"(after {res.attempts} attempt(s))")
        logger.error(f"-{len(failed)} of {len(tasks)} solver runs failed. Fix the cause and run "
                     "abics_mlref again: finished configurations are not rerun.")
        sys.exit(1)

    # --- per-replica result files, as in the staged mode ------------------
    by_replica: Dict[int, List[TaskResult]] = {r: [] for r in range(nreplicas)}
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
