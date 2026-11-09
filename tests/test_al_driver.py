import json
import os
import sys
import tempfile
import unittest
import unittest.mock

from abics.applications.latgas_abinitio_interface import al_driver as drv


class TestNextAction(unittest.TestCase):
    def test_sequence(self):
        # a fresh directory
        self.assertEqual(drv.next_action([], 0, 2), ("mlref", 0))
        # after each stage
        self.assertEqual(drv.next_action(["AL0"], 0, 2), ("train", 0))
        self.assertEqual(drv.next_action(["AL0", "train"], 0, 2), ("sampling", 0))
        self.assertEqual(drv.next_action(["AL0", "train", "MC0"], 1, 2), ("mlref", 1))
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1"], 1, 2), ("train", 1))
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1", "train"], 1, 2), ("sampling", 1))
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1", "train", "MC1"], 2, 2), ("done", 2))
        # raising max_iterations continues a finished loop
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1", "train", "MC1"], 2, 3), ("mlref", 2))
        # lowering it below the completed count is simply "done"
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1", "train", "MC1"], 2, 1), ("done", 2))
        # ... also in the middle of an iteration: the stages beyond the target are not continued
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1"], 1, 1), ("done", 1))
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1", "train"], 1, 1), ("done", 1))
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1", "train"], 2, 1), ("done", 1))  # interrupted MC1
        self.assertEqual(drv.completed_mc_runs(["AL0", "train", "MC0", "AL1", "train"]), 1)

    def test_interrupted_sampling_is_recognized(self):
        # abics_sampling creates MC<n> first and records it only when it finishes
        self.assertEqual(drv.next_action(["AL0", "train"], 1, 2), ("sampling_interrupted", 0))
        self.assertEqual(drv.next_action(["AL0", "train", "MC0", "AL1", "train"], 2, 3), ("sampling_interrupted", 1))

    def test_inconsistent_state_is_refused(self):
        with self.assertRaises(drv.DriverError):
            drv.next_action([], 1, 2)  # MC0 exists without a progress file
        with self.assertRaises(drv.DriverError):
            drv.next_action(["AL0"], 1, 2)  # MC0 exists but the progress says AL0 only
        with self.assertRaises(drv.DriverError):
            drv.next_action(["AL0", "train", "MC0"], 2, 3)  # MC1 exists but is not recorded
        with self.assertRaises(drv.DriverError):
            drv.next_action(["AL0", "train"], 2, 3)  # two MC dirs after a single AL0/train
        with self.assertRaises(drv.DriverError):
            drv.next_action(["train"], 0, 2)  # 'train' without any AL entry
        with self.assertRaises(drv.DriverError):
            drv.next_action(["something"], 0, 2)

    def test_set_aside_interrupted_mc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.makedirs(os.path.join(tmpdir, "MC1", "0"))
            driver_dir = os.path.join(tmpdir, drv.DRIVER_DIR)
            os.makedirs(driver_dir)
            moved = drv.set_aside_interrupted_mc(tmpdir, 1, driver_dir)
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "MC1")))
            self.assertTrue(os.path.isdir(os.path.join(moved, "0")))
            self.assertTrue(os.path.basename(moved).startswith("MC1.interrupted-"))


class TestParamsAndCommands(unittest.TestCase):
    def test_driver_params(self):
        p = drv.DriverParams.from_dict(
            {"max_iterations": 4, "walltime": 3600, "resubmit_command": "sbatch job.sh",
             "stage_estimates": {"mlref": 600, "sampling": 1200}}
        )
        self.assertEqual((p.max_iterations, p.walltime, p.resubmit_command), (4, 3600.0, "sbatch job.sh"))
        self.assertEqual(p.stage_estimates, {"mlref": 600.0, "sampling": 1200.0})
        self.assertIn("{input}", p.train_command)
        self.assertIn("{nprocs}", p.sampling_command)
        for bad in ({"max_iterations": 0}, {"walltime": -1}, {"stage_estimates": {"dft": 10}},
                    {"stage_estimates": {"train": -1}}):
            with self.assertRaises(ValueError):
                drv.DriverParams.from_dict(bad)

    def test_expand_command_keeps_other_braces_and_expands_paths(self):
        with unittest.mock.patch.dict(os.environ, {"ABICS_TEST_BIN": "/opt/bin"}):
            args = drv.expand_command(
                "srun -n {nprocs} --cpu-bind=map_cpu:{0,1} $ABICS_TEST_BIN/abics_sampling {input} ~/x",
                nprocs=8, input="in.toml",
            )
        self.assertEqual(args[:4], ["srun", "-n", "8", "--cpu-bind=map_cpu:{0,1}"])
        self.assertEqual(args[4], "/opt/bin/abics_sampling")
        self.assertEqual(args[5], "in.toml")
        self.assertEqual(args[6], os.path.expanduser("~/x"))
        # unknown placeholders are left alone instead of raising
        self.assertEqual(drv.expand_command("echo {unknown}", input="x"), ["echo", "{unknown}"])
        # a value with spaces stays one argument
        self.assertEqual(drv.expand_command("abics_train {input}", input="/a dir/in.toml"), ["abics_train", "/a dir/in.toml"])

    def test_sampling_nprocs(self):
        self.assertEqual(drv.sampling_nprocs({"sampling": {"nreplicas": 4, "solver": {"base_input_dir": "./b"}}}), (4, 4))
        self.assertEqual(
            drv.sampling_nprocs({"sampling": {"nreplicas": 4, "solver": {"base_input_dir": ["./b1", "./b2", "./b3"]}}}),
            (4, 12),
        )


class TestWallClock(unittest.TestCase):
    def test_wall_clock(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            params = drv.DriverParams.from_dict({"walltime": 100, "walltime_margin": 10})
            clock = drv.WallClock(params, tmpdir)
            ok, message = clock.can_start("mlref")  # unknown duration: started with a warning
            self.assertTrue(ok)
            self.assertIn("unknown", message)
            clock.record("mlref", 95.0)  # the last mlref took 95 s: it would not fit any more
            ok, message = clock.can_start("mlref")
            self.assertFalse(ok)
            self.assertIn("95 s", message)
            clock.record("mlref", 20.0)  # the longest run is kept, not the last one
            self.assertEqual(clock.expected("mlref"), 95.0)
            # durations persist across runs
            with open(os.path.join(tmpdir, drv.DURATIONS_FILE)) as f:
                self.assertEqual(json.load(f), {"mlref": 95.0})
            self.assertEqual(drv.WallClock(params, tmpdir).durations, {"mlref": 95.0})

    def test_estimates_are_used_until_the_stage_has_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            params = drv.DriverParams.from_dict(
                {"walltime": 100, "walltime_margin": 10, "stage_estimates": {"sampling": 200}}
            )
            clock = drv.WallClock(params, tmpdir)
            ok, message = clock.can_start("sampling")
            self.assertFalse(ok)
            self.assertIn("200 s", message)
            clock.record("sampling", 30.0)
            self.assertEqual(clock.expected("sampling"), 30.0)
            self.assertTrue(clock.can_start("sampling")[0])
            unlimited = drv.WallClock(drv.DriverParams.from_dict({}), tmpdir)
            self.assertEqual(unlimited.can_start("mlref"), (True, ""))


class TestStages(unittest.TestCase):
    def test_run_child(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log = os.path.join(tmpdir, "x.log")
            self.assertEqual(drv.run_child([sys.executable, "-c", "print(42)"], log, 0), 0)
            self.assertEqual(drv.run_child([sys.executable, "-c", "import sys; sys.exit(3)"], log, 0), 3)
            with open(log) as f:
                content = f.read()
            self.assertIn("42", content)
            self.assertEqual(content.count("==="), 2)  # both runs are kept (appended)
            self.assertEqual(drv.run_child([sys.executable, "-c", "import time; time.sleep(30)"], log, 1), -1)
            self.assertEqual(drv.run_child(["no-such-command-abics", "x"], log, 0), 127)

    def test_wait_for_progress(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self.assertFalse(drv.wait_for_progress(tmpdir, "train", 0.3))
            with open(os.path.join(tmpdir, "ALloop.progress"), "w") as f:
                f.write("AL0\ntrain\n")
            self.assertTrue(drv.wait_for_progress(tmpdir, "train", 0.3))
            self.assertFalse(drv.wait_for_progress(tmpdir, "MC0", 0.3))

    def test_ensure_train_input_generates_aenet_files(self):
        params_root = {
            "sampling": {"solver": {"type": "aenet", "ignore_species": ["O"]}},
            "train": {"type": "aenet", "ignore_species": ["O"]},
            "config": {
                "unitcell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "supercell": [1, 1, 1],
                "base_structure": [{"type": "O", "coords": [[0, 0, 0]]}],
                "defect_structure": [
                    {"coords": [[0.5, 0.5, 0.5]], "groups": [{"name": "Al", "num": 1}, {"name": "Mg", "num": 0}]}
                ],
            },
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            with unittest.mock.patch.dict(os.environ, {"ABICS_TEST_TMP": tmpdir}):
                base = os.path.join(tmpdir, "aenet_train_input")
                params_root["train"]["base_input_dir"] = "$ABICS_TEST_TMP/aenet_train_input"  # expanded
                drv.ensure_train_input(params_root)
                self.assertTrue(os.path.isfile(os.path.join(base, "train", "train.in")))
                # an existing, non-empty directory is left alone
                with open(os.path.join(base, "train", "train.in"), "w") as f:
                    f.write("custom\n")
                drv.ensure_train_input(params_root)
                with open(os.path.join(base, "train", "train.in")) as f:
                    self.assertEqual(f.read(), "custom\n")
            # other trainers are not touched
            params_root["train"]["type"] = "nequip"
            params_root["train"]["base_input_dir"] = os.path.join(tmpdir, "nequip_input")
            drv.ensure_train_input(params_root)
            self.assertFalse(os.path.exists(params_root["train"]["base_input_dir"]))


class TestSummaryAndMetadata(unittest.TestCase):
    def write_step(self, rootdir, n, rows_by_replica, stats):
        for r, rows in enumerate(rows_by_replica):
            os.makedirs(os.path.join(rootdir, f"AL{n}", str(r)))
            with open(os.path.join(rootdir, f"AL{n}", str(r), "energy_corr.dat"), "w") as f:
                for row in rows:
                    f.write(" ".join(str(x) for x in row) + "\n")
        if stats is not None:
            with open(os.path.join(rootdir, f"AL{n}", "convergence.dat"), "w") as f:
                f.write("# header\n# replica  n  rmse  mae  max_error  max_error_sample\n")
                for r in range(len(rows_by_replica)):
                    f.write(f"{r}  2  0.1  0.1  0.2  0\n")
                f.write(f"all  {stats[0]}  {stats[1]}  {stats[2]}  {stats[3]}  0\n")

    def test_summary_and_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self.write_step(tmpdir, 0, [[("NaN", -1.0, 0), ("NaN", -2.0, 1)]] * 2, None)
            self.write_step(tmpdir, 1, [[(-1.0, -1.1, 0), (-2.0, -2.2, 1)]] * 2, (4, 0.158, 0.15, 0.2))
            self.write_step(tmpdir, 2, [[(-1.0, -1.05, 0), (-2.0, -2.1, 1)]] * 2, (4, 0.079, 0.075, 0.1))
            self.write_step(tmpdir, 3, [[(-1.0, -1.05, 0)]] * 2, None)  # AL3 without convergence.dat (in progress)

            lines = drv.write_summary(tmpdir)
            data = [l.split() for l in lines if not l.startswith("#")]
            self.assertEqual(data, [["1", "4", "0.158", "0.15", "0.2"], ["2", "4", "0.079", "0.075", "0.1"]])
            self.assertTrue(os.path.isfile(os.path.join(tmpdir, drv.SUMMARY_FILE)))

            self.assertEqual(drv.count_reference_data(tmpdir, 3), 12)  # AL0..AL2 only
            self.assertEqual(drv.count_reference_data(tmpdir, 4), 14)

            tomlfile = os.path.join(tmpdir, "input.toml")
            with open(tomlfile, "w") as f:
                f.write("[al]\nmax_iterations = 3\n")
            params_root = {
                "sampling": {"nreplicas": 2, "kTstart": 1000.0, "kTend": 1500.0,
                             "solver": {"type": "aenet", "base_input_dir": "./baseinput", "ignore_species": ["O"]}},
                "mlref": {"seed": 7, "solver": {"type": "qe", "base_input_dir": "./baseinput_ref", "perturb": 0.05}},
                "train": {"type": "aenet", "base_input_dir": "./aenet_train_input", "ignore_species": ["O"]},
                "config": {"supercell": [1, 1, 1],
                           "base_structure": [{"type": "O"}],
                           "defect_structure": [{"groups": [{"name": "Al"}, {"name": "Mg"}]}]},
            }
            driver_dir = os.path.join(tmpdir, drv.DRIVER_DIR)
            os.makedirs(driver_dir)
            path = drv.write_model_metadata(tmpdir, params_root, tomlfile, 3, driver_dir)
            with open(path) as f:
                meta = json.load(f)
            self.assertEqual(meta["provenance"]["iterations"], 3)
            self.assertEqual(meta["provenance"]["reference_configurations"], 12)
            self.assertEqual(meta["provenance"]["reference_solver"]["type"], "qe")
            self.assertEqual(meta["provenance"]["mlref_seed"], 7)
            self.assertIsNone(meta["provenance"]["base_model"])
            self.assertEqual(meta["scope"]["species"], ["Al", "Mg", "O"])
            self.assertEqual(meta["scope"]["kT_range"], [1000.0, 1500.0])
            self.assertEqual(meta["model_files"], "./baseinput")
            self.assertEqual(sorted(meta["accuracy"]), ["AL1", "AL2"])
            self.assertAlmostEqual(meta["accuracy"]["AL2"]["rmse"], 0.079)
            self.assertEqual(len(meta["input_sha256"]), 64)


class TestRunDriverWithStubbedStages(unittest.TestCase):
    """The loop itself, with the three stages replaced by functions that write what the programs write."""

    params_root = {
        "sampling": {"nreplicas": 2, "solver": {"type": "aenet", "base_input_dir": "./baseinput"}},
        "mlref": {"nreplicas": 2, "ndata": 2, "solver": {"type": "user", "path": "./mock.py",
                                                         "run_scheme": "subprocess", "function": "m.f"}},
        "train": {"type": "aenet", "base_input_dir": "./aenet_train_input"},
        "config": {"supercell": [1, 1, 1], "base_structure": [], "defect_structure": []},
    }

    def stubs(self, tmpdir, calls, fail_sampling_at=None, slow=None):
        def mlref(params_root):  # stands in for run_integrated_step
            n = drv.newest_mc_index(tmpdir)
            calls.append(("mlref", n))
            self.assertEqual(params_root["mlref"]["exec_mode"], "integrated")
            os.makedirs(os.path.join(tmpdir, f"AL{n}", "0"), exist_ok=True)
            al_path = os.path.join(tmpdir, al_progress)
            with open(al_path, "a") as f:
                f.write(f"AL{n}\n")
            return n

        def train(params, params_root, tomlfile, n, driver_dir, rootdir):
            calls.append(("train", n))
            with open(os.path.join(tmpdir, al_progress), "a") as f:
                f.write("train\n")

        def sampling(params, params_root, tomlfile, n, driver_dir, rootdir):
            calls.append(("sampling", n))
            os.makedirs(os.path.join(tmpdir, f"MC{n}"), exist_ok=True)  # created first, as abics_sampling does
            if fail_sampling_at == n:
                raise drv.DriverError(f"MC{n}: sampling failed (test)")
            with open(os.path.join(tmpdir, al_progress), "a") as f:
                f.write(f"MC{n}\n")

        al_progress = "ALloop.progress"
        return mlref, train, sampling

    def run_loop(self, tmpdir, max_iterations, calls, **kwargs):
        tomlfile = os.path.join(tmpdir, "input.toml")
        if not os.path.exists(tomlfile):
            with open(tomlfile, "w") as f:
                f.write("[al]\n")
        mlref, train, sampling = self.stubs(tmpdir, calls, **kwargs)
        cwd = os.getcwd()
        os.chdir(tmpdir)
        try:
            with unittest.mock.patch.object(drv, "run_integrated_step", mlref), \
                 unittest.mock.patch.object(drv, "run_train", train), \
                 unittest.mock.patch.object(drv, "run_sampling", sampling), \
                 unittest.mock.patch.object(drv, "require_libensemble", lambda: None):
                return drv.run_driver(self.params_root, tomlfile, max_iterations)
        finally:
            os.chdir(cwd)

    def test_full_loop_then_extension(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            calls = []
            self.assertEqual(self.run_loop(tmpdir, 2, calls), "done")
            self.assertEqual(calls, [("mlref", 0), ("train", 0), ("sampling", 0),
                                     ("mlref", 1), ("train", 1), ("sampling", 1)])
            self.assertEqual(drv.read_progress(tmpdir), ["AL0", "train", "MC0", "AL1", "train", "MC1"])
            with open(os.path.join(tmpdir, drv.DRIVER_DIR, drv.METADATA_FILE)) as f:
                self.assertEqual(json.load(f)["provenance"]["iterations"], 2)
            # extension: only the third iteration runs
            calls.clear()
            self.assertEqual(self.run_loop(tmpdir, 3, calls), "done")
            self.assertEqual(calls, [("mlref", 2), ("train", 2), ("sampling", 2)])
            # a lower target than what exists: nothing runs, metadata reports what exists
            calls.clear()
            self.assertEqual(self.run_loop(tmpdir, 1, calls), "done")
            self.assertEqual(calls, [])
            with open(os.path.join(tmpdir, drv.DRIVER_DIR, drv.METADATA_FILE)) as f:
                self.assertEqual(json.load(f)["provenance"]["iterations"], 3)

    def test_interrupted_sampling_is_moved_aside_and_repeated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            calls = []
            with self.assertRaises(drv.DriverError):
                self.run_loop(tmpdir, 2, calls, fail_sampling_at=0)
            # MC0 exists, progress ends with "train": the state of a killed sampling run
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "MC0")))
            self.assertEqual(drv.read_progress(tmpdir)[-1], "train")
            calls.clear()
            self.assertEqual(self.run_loop(tmpdir, 1, calls), "done")
            self.assertEqual(calls, [("sampling", 0)])  # repeated, nothing else
            moved = [d for d in os.listdir(os.path.join(tmpdir, drv.DRIVER_DIR)) if d.startswith("MC0.interrupted-")]
            self.assertEqual(len(moved), 1)

    def run_with_al(self, tmpdir, al_section, max_iterations=2, **stub_kwargs):
        tomlfile = os.path.join(tmpdir, "input.toml")
        with open(tomlfile, "w") as f:
            f.write("[al]\n")
        params_root = dict(self.params_root)
        params_root["al"] = al_section
        calls = []
        mlref, train, sampling = self.stubs(tmpdir, calls, **stub_kwargs)
        cwd = os.getcwd()
        os.chdir(tmpdir)
        try:
            with unittest.mock.patch.object(drv, "run_integrated_step", mlref), \
                 unittest.mock.patch.object(drv, "run_train", train), \
                 unittest.mock.patch.object(drv, "run_sampling", sampling), \
                 unittest.mock.patch.object(drv, "require_libensemble", lambda: None):
                return drv.run_driver(params_root, tomlfile, max_iterations), calls
        finally:
            os.chdir(cwd)

    def test_walltime_stop_runs_resubmit_command(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            marker = os.path.join(tmpdir, "resubmitted")
            al = {"walltime": 100, "walltime_margin": 10, "stage_estimates": {"train": 1000},
                  "resubmit_command": f"{sys.executable} -c \"open('{marker}', 'w').close()\""}
            state, calls = self.run_with_al(tmpdir, al)
            # mlref (unknown duration) ran; the training (estimated 1000 s) does not fit into 100 s
            self.assertEqual(state, "walltime")
            self.assertEqual(calls, [("mlref", 0)])
            self.assertTrue(os.path.exists(marker))

    def test_no_progress_means_error_instead_of_resubmission(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            marker = os.path.join(tmpdir, "resubmitted")
            al = {"walltime": 100, "walltime_margin": 10, "stage_estimates": {"mlref": 1000},
                  "resubmit_command": f"{sys.executable} -c \"open('{marker}', 'w').close()\""}
            with self.assertRaises(drv.DriverError) as cm:
                self.run_with_al(tmpdir, al)
            self.assertIn("does not fit into the wall time even at the start", str(cm.exception))
            self.assertFalse(os.path.exists(marker))  # no endless resubmission loop

    def test_failing_resubmit_command_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            al = {"walltime": 100, "walltime_margin": 10, "stage_estimates": {"train": 1000},
                  "resubmit_command": f"{sys.executable} -c 'import sys; sys.exit(4)'"}
            with self.assertRaises(drv.DriverError) as cm:
                self.run_with_al(tmpdir, al)
            self.assertIn("exited with status 4", str(cm.exception))

    def test_reload_is_refused(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            params_root = {k: dict(v) for k, v in self.params_root.items()}
            params_root["sampling"]["reload"] = True
            with unittest.mock.patch.object(drv, "require_libensemble", lambda: None):
                with self.assertRaises(drv.DriverError) as cm:
                    drv.run_driver(params_root, os.path.join(tmpdir, "input.toml"), 1)
            self.assertIn("reload", str(cm.exception))

    def test_lock_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            driver_dir = os.path.join(tmpdir, drv.DRIVER_DIR)
            os.makedirs(driver_dir)
            # a live process (this one) holds the lock
            with open(os.path.join(driver_dir, drv.LOCK_FILE), "w") as f:
                f.write(f"{os.uname().nodename} {os.getpid()} now\n")
            with self.assertRaises(drv.DriverError):
                drv.DriverLock(driver_dir).acquire()
            # a stale lock (dead pid) is replaced
            with open(os.path.join(driver_dir, drv.LOCK_FILE), "w") as f:
                f.write(f"{os.uname().nodename} 999999999 then\n")
            lock = drv.DriverLock(driver_dir)
            lock.acquire()
            with open(lock.path) as f:
                self.assertIn(str(os.getpid()), f.read())
            lock.release()
            self.assertFalse(os.path.exists(lock.path))
            # the driver removes its lock even when a stage fails
            with self.assertRaises(drv.DriverError):
                self.run_with_al(tmpdir, {}, fail_sampling_at=0)
            self.assertFalse(os.path.exists(lock.path))


if __name__ == "__main__":
    unittest.main()
