import os
import tempfile
import unittest
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
            with open(os.path.join(task.input_dir, "baseinput0", "energy.dat"), "w") as f:
                f.write("-3.0\n")

            def launch(run_dir, attempt):
                raise AssertionError("must not be called")

            status, energy, _, attempts = ai.run_task_stages(task, FakeSolver(), alparams(), 1, launch)
        self.assertEqual((status, energy, attempts), ("ok", -3.0, 0))

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
        with self.assertRaises(ValueError):
            ai.IntegratedParams.from_dict({"nworkers": 0})

    def test_copy_tree_filtered(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            src = os.path.join(tmpdir, "src")
            os.makedirs(os.path.join(src, "sub"))
            for name in ("a.out", "b.wfc1", "sub/c.dat", "sub/d.wfc2"):
                with open(os.path.join(src, name), "w") as f:
                    f.write(name)
            dst = os.path.join(tmpdir, "dst")
            ai.copy_tree_filtered(src, dst, ["*.wfc*"])
            copied = sorted(
                os.path.relpath(os.path.join(r, f), dst) for r, _, fs in os.walk(dst) for f in fs
            )
            self.assertEqual(copied, ["a.out", "sub/c.dat"])

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
