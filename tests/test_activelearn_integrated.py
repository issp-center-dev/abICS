import importlib.util
import os
import tempfile
import time
import unittest
import unittest.mock
from collections import namedtuple

import numpy as np
from pymatgen.core import Lattice, Structure

from abics.applications.latgas_abinitio_interface import activelearn as al
from abics.applications.latgas_abinitio_interface import activelearn_integrated as ai


class FakeSolver:
    """A solver whose 'run' is the launcher writing energy.dat; parsing reads it back."""

    Phys = namedtuple("PhysValues", ("energy", "structure"))

    class Input:
        def __init__(self):
            self.st = None
            self.base_dirs = []

        def from_directory(self, base_input_dir):
            self.base_dirs.append(base_input_dir)

        def update_info_by_structure(self, structure):
            self.st = structure

        def write_input(self, output_dir):
            os.makedirs(output_dir, exist_ok=True)
            self.st.to(fmt="POSCAR", filename=os.path.join(output_dir, "structure.vasp"))

        def cl_args(self, nprocs, nthreads, workdir):
            return [workdir]

    class Output:
        def get_results(self, workdir):
            with open(os.path.join(workdir, "energy.dat")) as f:
                ene = float(f.read())
            st = Structure.from_file(os.path.join(workdir, "structure.vasp"))
            return FakeSolver.Phys(ene, st)

    def __init__(self):
        self.input = FakeSolver.Input()
        self.output = FakeSolver.Output()


def make_task(tmpdir, replica=0, index=3):
    lattice = Lattice.cubic(4.0)
    st = Structure(lattice, ["Al", "Mg"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    input_dir = os.path.join(tmpdir, f"AL1/{replica}/input{index:02d}")
    stage0 = os.path.join(input_dir, "baseinput0")
    os.makedirs(stage0)
    st.to(fmt="POSCAR", filename=os.path.join(stage0, "structure.vasp"))
    st.to(fmt="POSCAR", filename=os.path.join(input_dir, "initial.vasp"))
    return ai.ReferenceTask(replica, index, input_dir)


def alparams(nstages=1, ignore_species=None):
    params = ai.ALParams()
    params.base_input_dir = ["./base"] * nstages
    params.ignore_species = ignore_species
    return params


class TestRunTaskStages(unittest.TestCase):
    def test_single_stage_ok(self):
        calls = []

        def launch(run_dir, attempt):
            calls.append((os.path.basename(run_dir), attempt))
            with open(os.path.join(run_dir, "energy.dat"), "w") as f:
                f.write("-1.5\n")
            return "ok"

        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir)
            status, energy, st, attempts = ai.run_task_stages(task, FakeSolver(), alparams(), 1, launch)
        self.assertEqual(status, "ok")
        self.assertEqual(energy, -1.5)
        self.assertEqual(len(st), 2)
        self.assertEqual(attempts, 1)
        self.assertEqual(calls, [("baseinput0", 1)])

    def test_retry_then_failure(self):
        def launch(run_dir, attempt):
            return f"solver exited with status {attempt}"

        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir)
            status, energy, st, attempts = ai.run_task_stages(task, FakeSolver(), alparams(), 2, launch)
        self.assertEqual(status, "stage 0: solver exited with status 3")
        self.assertTrue(np.isnan(energy))
        self.assertIsNone(st)
        self.assertEqual(attempts, 3)

    def test_ok_launch_but_unparseable_output_is_retried(self):
        attempts_seen = []

        def launch(run_dir, attempt):
            attempts_seen.append(attempt)
            if attempt == 2:  # the first run "finished" without writing energy.dat
                with open(os.path.join(run_dir, "energy.dat"), "w") as f:
                    f.write("-2.0\n")
            return "ok"

        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir)
            status, energy, _, attempts = ai.run_task_stages(task, FakeSolver(), alparams(), 1, launch)
        self.assertEqual(status, "ok")
        self.assertEqual(energy, -2.0)
        self.assertEqual(attempts_seen, [1, 2])
        self.assertEqual(attempts, 2)

    def test_finished_stage_is_not_rerun(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir)
            stage0 = os.path.join(task.input_dir, "baseinput0")
            with open(os.path.join(stage0, "energy.dat"), "w") as f:
                f.write("-3.0\n")
            ai.mark_stage_done(stage0, 1)

            def launch(run_dir, attempt):
                raise AssertionError("must not be called")

            status, energy, _, attempts = ai.run_task_stages(task, FakeSolver(), alparams(), 1, launch)
        self.assertEqual((status, energy, attempts), ("ok", -3.0, 0))

    def test_output_without_marker_is_rerun(self):
        # output of an interrupted run: parseable but not marked complete
        calls = []

        def launch(run_dir, attempt):
            calls.append(attempt)
            with open(os.path.join(run_dir, "energy.dat"), "w") as f:
                f.write("-4.0\n")
            return "ok"

        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir)
            stage0 = os.path.join(task.input_dir, "baseinput0")
            with open(os.path.join(stage0, "energy.dat"), "w") as f:
                f.write("-3.0\n")  # stale partial output
            status, energy, _, attempts = ai.run_task_stages(task, FakeSolver(), alparams(), 1, launch)
            self.assertTrue(ai.stage_done(stage0))
        self.assertEqual((status, energy, attempts, calls), ("ok", -4.0, 1, [1]))

    def test_marker_is_written_only_after_success(self):
        def launch(run_dir, attempt):
            return "solver exited with status 1"

        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir)
            status, _, _, _ = ai.run_task_stages(task, FakeSolver(), alparams(), 0, launch)
            self.assertFalse(ai.stage_done(os.path.join(task.input_dir, "baseinput0")))
        self.assertTrue(status.startswith("stage 0: solver exited"))

    def test_output_appearing_late_is_waited_for(self):
        import threading

        def launch(run_dir, attempt):
            def write_later():
                time.sleep(0.5)
                with open(os.path.join(run_dir, "energy.dat"), "w") as f:
                    f.write("-5.0\n")

            threading.Thread(target=write_later).start()
            return "ok"

        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir)
            status, energy, _, attempts = ai.run_task_stages(
                task, FakeSolver(), alparams(), 0, launch, output_wait=5.0
            )
        self.assertEqual((status, energy, attempts), ("ok", -5.0, 1))

    def test_parser_calling_exit_is_treated_as_not_finished(self):
        class ExitingOutput:
            def get_results(self, workdir):
                exit(1)  # as the QE/aenet/MLIP-3 parsers do on a missing output

        with tempfile.TemporaryDirectory() as tmpdir:
            self.assertIsNone(ai.try_parse(ExitingOutput(), tmpdir))

    def test_two_stages_are_chained(self):
        stages = []

        def launch(run_dir, attempt):
            stages.append(os.path.basename(run_dir))
            with open(os.path.join(run_dir, "energy.dat"), "w") as f:
                f.write(f"-{len(stages)}.0\n")
            return "ok"

        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir)
            solver = FakeSolver()
            status, energy, _, attempts = ai.run_task_stages(task, solver, alparams(2), 0, launch)
            self.assertTrue(os.path.isfile(os.path.join(task.input_dir, "baseinput1", "structure.vasp")))
        self.assertEqual(status, "ok")
        self.assertEqual(energy, -2.0)
        self.assertEqual(stages, ["baseinput0", "baseinput1"])
        self.assertEqual(solver.input.base_dirs, ["./base"])  # the second stage's base input was read
        self.assertEqual(attempts, 2)


class TestFinishTask(unittest.TestCase):
    def test_structure_file_and_displacement(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir, replica=1, index=7)
            replica_dir = os.path.join(tmpdir, "AL1", "1")
            # the solver moved the second atom by 0.4 A along x (cubic 4 A cell)
            st = Structure.from_file(os.path.join(task.input_dir, "baseinput0", "structure.vasp"))
            st.translate_sites([1], [0.1, 0, 0], frac_coords=True)
            dmax_id, dmax = ai.finish_task(task, st, alparams(), 1, replica_dir)
            self.assertEqual(dmax_id, 1)
            self.assertAlmostEqual(dmax, 0.4)
            written = Structure.from_file(os.path.join(replica_dir, "structure.7.vasp"))
            self.assertEqual(len(written), 2)

    def test_ignored_species_are_restored_and_al0_sorted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            task = make_task(tmpdir, replica=0, index=2)
            replica_dir = os.path.join(tmpdir, "AL1", "0")
            # initial.vasp contains O, the solver's structure does not
            lattice = Lattice.cubic(4.0)
            full = Structure(lattice, ["Mg", "Al", "O"], [[0.5, 0.5, 0.5], [0, 0, 0], [0.25, 0.25, 0.25]])
            full.to(fmt="POSCAR", filename=os.path.join(task.input_dir, "initial.vasp"))
            st = Structure(lattice, ["Mg", "Al"], [[0.5, 0.5, 0.5], [0, 0, 0]])
            ai.finish_task(task, st, alparams(ignore_species=["O"]), 0, replica_dir)
            written = Structure.from_file(os.path.join(replica_dir, "structure.2.vasp"))
            self.assertEqual([s.species_string for s in written], ["Al", "Mg", "O"])


class TestHelpers(unittest.TestCase):
    def test_integrated_params(self):
        p = ai.IntegratedParams.from_dict({"nworkers": 4, "nprocs": 2, "timeout": 30, "scratch_dir": "/tmp/x"})
        self.assertEqual((p.nworkers, p.nprocs, p.timeout, p.retries, p.scratch_dir), (4, 2, 30.0, 1, "/tmp/x"))
        self.assertEqual(ai.IntegratedParams.from_dict({}).nworkers, 1)
        for bad in (
            {"nworkers": 0}, {"nthreads": 0}, {"timeout": -1}, {"retries": -1},
            {"nprocs": 4, "cores": 2},  # no run could be scheduled
            {"scratch_dir": "/scratch/with space"},
        ):
            with self.assertRaises(ValueError):
                ai.IntegratedParams.from_dict(bad)
        self.assertEqual(ai.IntegratedParams.from_dict({"nprocs": 4, "cores": 8}).cores, 8)

    def test_scratch_dir_is_expanded(self):
        with unittest.mock.patch.dict(os.environ, {"ABICS_TEST_SCRATCH": "/scr"}):
            p = ai.IntegratedParams.from_dict({"scratch_dir": "$ABICS_TEST_SCRATCH/$USER"})
        self.assertTrue(p.scratch_dir.startswith("/scr/"))
        self.assertNotIn("$", p.scratch_dir)
        self.assertEqual(ai.IntegratedParams.from_dict({"scratch_dir": "~"}).scratch_dir, os.path.expanduser("~"))

    def test_launched_by_mpi(self):
        clean = {k: v for k, v in os.environ.items() if k not in ai.MPI_LAUNCHER_VARIABLES}
        with unittest.mock.patch.dict(os.environ, clean, clear=True):
            self.assertFalse(ai.launched_by_mpi())
            # even "mpirun -np 1" is refused: the nested mpirun of the solver would fail
            with unittest.mock.patch.dict(os.environ, {"OMPI_COMM_WORLD_SIZE": "1"}):
                self.assertTrue(ai.launched_by_mpi())
            with unittest.mock.patch.dict(os.environ, {"PMI_SIZE": "4"}):
                self.assertTrue(ai.launched_by_mpi())

    def test_seed_zero_gives_independent_replicas(self):
        ai.seed_replica(0, 0)
        a = np.random.random(4)
        ai.seed_replica(0, 1)
        b = np.random.random(4)
        self.assertFalse(np.array_equal(a, b))

    def test_check_resumable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            replica_dirs = [os.path.join(tmpdir, "AL1", "0"), os.path.join(tmpdir, "AL1", "1")]
            fmt = [al.input_dir_format(80)] * 2
            for name in ("input00", "input08", "input17"):
                os.makedirs(os.path.join(replica_dirs[0], name))
            # replica 1 was never generated
            energy_ref = [{0: -1.0, 8: -2.0, 17: -3.0}, {}]
            tasks, problems = ai.check_resumable(replica_dirs, fmt, 3, energy_ref)
            self.assertEqual([(t.replica, t.index) for t in tasks], [(0, 0), (0, 8), (0, 17)])
            self.assertEqual(len(problems), 1)
            self.assertIn("does not exist", problems[0])

            os.makedirs(os.path.join(replica_dirs[1], "input00"))
            tasks, problems = ai.check_resumable(replica_dirs, fmt, 3, energy_ref)
            self.assertEqual(len(problems), 2)  # wrong count, and index 0 not in MC samples of replica 1
            self.assertTrue(any("expected ndata = 3" in p for p in problems))
            self.assertTrue(any("no sample" in p for p in problems))

            # a complete AL0 (no energy_ref) resumes without problems
            tasks, problems = ai.check_resumable(replica_dirs[:1], fmt[:1], 3, None)
            self.assertEqual(problems, [])
            self.assertEqual(len(tasks), 3)

            # the directory names must match the current zero-padding width
            tasks, problems = ai.check_resumable(replica_dirs[:1], [al.input_dir_format(1000)], 3, None)
            self.assertEqual(tasks, [])
            self.assertEqual(len(problems), 3)
            self.assertIn("is not named input0000", problems[0])


@unittest.skipUnless(importlib.util.find_spec("libensemble"), "libEnsemble not installed")
class TestRunReferenceTask(unittest.TestCase):
    """The worker-side function, with a fake libEnsemble executor."""

    class FakeTask:
        def __init__(self, run_dir, succeed):
            self.run_dir = run_dir
            self.succeed = succeed
            self.success = False
            self.errcode = None

        def wait(self, timeout=None):
            # "the solver": computes the energy in the current directory (= run_dir)
            assert os.path.realpath(os.getcwd()) == os.path.realpath(self.run_dir)
            if self.succeed:
                with open(os.path.join(self.run_dir, "energy.dat"), "w") as f:
                    f.write("-7.5\n")
                open(os.path.join(self.run_dir, "machinefile_autogen_for_worker_1_task_0"), "w").close()
            self.success = self.succeed
            self.errcode = 0 if self.succeed else 5

        def kill(self):
            pass

    class FakeExecutor:
        def __init__(self, succeed=True):
            self.succeed = succeed
            self.submits = []

        def submit(self, **kwargs):
            self.submits.append(kwargs)
            return TestRunReferenceTask.FakeTask(os.getcwd(), self.succeed)

    def run_worker(self, tmpdir, fake, scratch_root="", al_step=1):
        from libensemble.executors.executor import Executor

        task = make_task(tmpdir, replica=0, index=3)
        user = {
            "tasks": [{"replica": 0, "index": 3, "input_dir": task.input_dir}],
            "mlref_solver": {"type": "user", "function": "fakemock.calc_energy",
                             "path": "/bin/true", "run_scheme": "subprocess", "perturb": 0.0},
            "integrated": {"nworkers": 1, "nprocs": 1, "nthreads": 1, "timeout": 0.0, "retries": 0,
                           "output_wait": 0.0, "scratch_dir": scratch_root and os.path.dirname(scratch_root),
                           "copy_back_exclude": [], "cores": 0},
            "rootdir": tmpdir,
            "replica_dirs": [os.path.join(tmpdir, "AL1", "0")],
            "al_step": al_step,
            "scratch_root": scratch_root,
        }
        # the user-function solver imports its module from rootdir
        with open(os.path.join(tmpdir, "fakemock.py"), "w") as f:
            f.write("def calc_energy(st):\n    return 0.0\n")
        H = np.zeros(1, dtype=[("task_index", int)])
        sim_specs = {
            "out": [("task_index", int), ("status", "U256"), ("energy", float),
                    ("dmax_id", int), ("dmax", float), ("attempts", int)],
            "user": user,
        }
        cwd = os.getcwd()
        try:
            with unittest.mock.patch.object(Executor, "executor", fake, create=True):
                H_o, _ = ai.run_reference_task(H, {}, sim_specs, {})
        finally:
            os.chdir(cwd)
        return task, H_o

    def test_success_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fake = self.FakeExecutor(succeed=True)
            task, H_o = self.run_worker(tmpdir, fake)
            stage0 = os.path.join(task.input_dir, "baseinput0")
            self.assertEqual(str(H_o["status"][0]), "ok")
            self.assertEqual(float(H_o["energy"][0]), -7.5)
            self.assertEqual(int(H_o["attempts"][0]), 1)
            self.assertTrue(ai.stage_done(stage0))
            self.assertTrue(os.path.isfile(os.path.join(tmpdir, "AL1", "0", "structure.3.vasp")))
            self.assertEqual(fake.submits[0]["app_args"], stage0)  # cl_args of the user solver
            self.assertFalse(any(n.startswith("machinefile") for n in os.listdir(stage0)))

    def test_solver_failure_is_reported(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            task, H_o = self.run_worker(tmpdir, self.FakeExecutor(succeed=False))
            self.assertEqual(str(H_o["status"][0]), "stage 0: solver exited with status 5")
            self.assertTrue(np.isnan(H_o["energy"][0]))
            self.assertFalse(ai.stage_done(os.path.join(task.input_dir, "baseinput0")))

    def test_internal_error_becomes_a_status(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fake = self.FakeExecutor(succeed=True)
            task_dir = os.path.join(tmpdir, "AL1", "0", "input03")
            os.makedirs(os.path.join(task_dir, "baseinput0"))
            # make_task would write initial.vasp; here it is missing -> finish_task fails
            with open(os.path.join(task_dir, "baseinput0", "structure.vasp"), "w") as f:
                f.write(Structure(Lattice.cubic(4.0), ["Al"], [[0, 0, 0]]).to(fmt="POSCAR"))
            os.makedirs(os.path.join(tmpdir, "AL1", "0"), exist_ok=True)
            from libensemble.executors.executor import Executor
            with open(os.path.join(tmpdir, "fakemock.py"), "w") as f:
                f.write("def calc_energy(st):\n    return 0.0\n")
            sim_specs = {
                "out": [("task_index", int), ("status", "U256"), ("energy", float),
                        ("dmax_id", int), ("dmax", float), ("attempts", int)],
                "user": {
                    "tasks": [{"replica": 0, "index": 3, "input_dir": task_dir}],
                    "mlref_solver": {"type": "user", "function": "fakemock.calc_energy",
                                     "path": "/bin/true", "run_scheme": "subprocess", "perturb": 0.0},
                    "integrated": {"nworkers": 1, "nprocs": 1, "nthreads": 1, "timeout": 0.0, "retries": 0,
                                   "output_wait": 0.0, "scratch_dir": "", "copy_back_exclude": [], "cores": 0},
                    "rootdir": tmpdir, "replica_dirs": [os.path.join(tmpdir, "AL1", "0")],
                    "al_step": 1, "scratch_root": "",
                },
            }
            cwd = os.getcwd()
            try:
                with unittest.mock.patch.object(Executor, "executor", fake, create=True):
                    H_o, _ = ai.run_reference_task(np.zeros(1, dtype=[("task_index", int)]), {}, sim_specs, {})
            finally:
                os.chdir(cwd)
            self.assertTrue(str(H_o["status"][0]).startswith("internal error: "))
            self.assertIn("initial.vasp", str(H_o["status"][0]))

    def test_scratch_copy_back_and_cleanup(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scratch_root = os.path.join(tmpdir, "scratch", "abics_mlref_AL1")
            os.makedirs(scratch_root)
            fake = self.FakeExecutor(succeed=True)
            task, H_o = self.run_worker(tmpdir, fake, scratch_root=scratch_root)
            stage0 = os.path.join(task.input_dir, "baseinput0")
            self.assertEqual(str(H_o["status"][0]), "ok")
            self.assertTrue(fake.submits[0]["stdout"].startswith(scratch_root))  # ran in scratch
            self.assertTrue(os.path.isfile(os.path.join(stage0, "energy.dat")))  # copied back
            self.assertFalse(any(n.startswith("machinefile") for n in os.listdir(stage0)))
            self.assertEqual(os.listdir(scratch_root), [])  # run directory removed

    def test_copy_tree_filtered(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            src = os.path.join(tmpdir, "src")
            os.makedirs(os.path.join(src, "sub"))
            os.makedirs(os.path.join(src, "pwscf.save", "deep"))
            for name in (
                "a.out", "b.wfc1", "sub/c.dat", "sub/d.wfc2", "machinefile_autogen_for_worker_1_task_0",
                "pwscf.save/charge-density.dat", "pwscf.save/deep/x",
            ):
                with open(os.path.join(src, name), "w") as f:
                    f.write(name)
            dst = os.path.join(tmpdir, "dst")
            ai.copy_tree_filtered(src, dst, ["*.wfc*", "*.save"])
            copied = sorted(
                os.path.relpath(os.path.join(r, f), dst) for r, _, fs in os.walk(dst) for f in fs
            )
            # "*.save" excludes the whole directory; launcher machinefiles are never copied
            self.assertEqual(copied, ["a.out", "sub/c.dat"])
            self.assertFalse(os.path.exists(os.path.join(dst, "pwscf.save")))

    def test_remove_launcher_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for name in ("machinefile_autogen_for_worker_2_task_5", "stdout"):
                open(os.path.join(tmpdir, name), "w").close()
            ai.remove_launcher_files(tmpdir)
            self.assertEqual(os.listdir(tmpdir), ["stdout"])

    def test_resolve_solver_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            exe = os.path.join(tmpdir, "solver.x")
            open(exe, "w").close()
            self.assertEqual(ai.resolve_solver_path(exe), exe)
            self.assertTrue(os.path.isabs(ai.resolve_solver_path("python3")))
            with self.assertRaises(FileNotFoundError):
                ai.resolve_solver_path("no-such-command-abics")

    def test_seed_replica_matches_rx_mpi_init(self):
        ai.seed_replica(100, 2)
        a = np.random.random()
        np.random.seed(100 + 2 * 137)
        self.assertEqual(a, np.random.random())

    def test_task_list_gen(self):
        gen_specs = {"out": [("task_index", int), ("resource_sets", int)], "user": {"ntasks": 3, "nprocs": 2}}
        H_o, _ = ai.task_list_gen(None, {}, gen_specs, {})
        self.assertEqual(list(H_o["task_index"]), [0, 1, 2])
        self.assertEqual(list(H_o["resource_sets"]), [2, 2, 2])


if __name__ == "__main__":
    unittest.main()
