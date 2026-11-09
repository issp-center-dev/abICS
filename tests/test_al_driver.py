import json
import os
import subprocess
import sys
import tempfile
import unittest
import unittest.mock

from abics.applications.latgas_abinitio_interface import al_driver as drv


def make_old(path, age=3600):
    """Set the modification times under ``path`` (two levels) ``age`` seconds into the past."""
    import time

    t = time.time() - age
    os.utime(path, (t, t))
    for entry in os.scandir(path):
        os.utime(entry.path, (t, t))
        if entry.is_dir():
            for sub in os.scandir(entry.path):
                os.utime(sub.path, (t, t))


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
            open(os.path.join(tmpdir, "MC1", "0", "obs.dat"), "w").close()
            driver_dir = os.path.join(tmpdir, drv.DRIVER_DIR)
            os.makedirs(driver_dir)
            # just written: a sampling run may still be alive
            with self.assertRaises(drv.DriverError) as cm:
                drv.set_aside_interrupted_mc(tmpdir, 1, driver_dir)
            self.assertIn("may still be active", str(cm.exception))
            make_old(os.path.join(tmpdir, "MC1"))
            moved = drv.set_aside_interrupted_mc(tmpdir, 1, driver_dir)
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "MC1")))
            self.assertTrue(os.path.isdir(os.path.join(moved, "0")))
            self.assertTrue(os.path.basename(moved).startswith("MC1.interrupted-"))

    def test_trained_steps(self):
        self.assertEqual(drv.trained_steps([]), 0)
        self.assertEqual(drv.trained_steps(["AL0"]), 0)
        self.assertEqual(drv.trained_steps(["AL0", "train"]), 1)
        self.assertEqual(drv.trained_steps(["AL0", "train", "MC0", "AL1"]), 1)
        self.assertEqual(drv.trained_steps(["AL0", "train", "MC0", "AL1", "train"]), 2)
        self.assertEqual(drv.trained_steps(["AL0", "train", "MC0", "AL1", "train", "MC1"]), 2)


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
        # a value with spaces stays one argument, quoted or not, and values are not expanded
        self.assertEqual(drv.expand_command("abics_train {input}", input="/a dir/in.toml"), ["abics_train", "/a dir/in.toml"])
        self.assertEqual(drv.expand_command('abics_train "{input}"', input="/a dir/in.toml"), ["abics_train", "/a dir/in.toml"])
        self.assertEqual(drv.expand_command("x {input}", input="$HOME/~")[1], "$HOME/~")

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

    def test_child_is_stopped_when_the_driver_is_interrupted(self):
        import signal

        class FakeProc:
            pid = 424242

            def __init__(self):
                self.waits = 0

            def wait(self, timeout=None):
                self.waits += 1
                if self.waits == 1:
                    raise KeyboardInterrupt
                return 0

            def poll(self):
                return 0 if self.waits >= 2 else None  # reaped by _stop_child's wait

        fake = FakeProc()
        signals = []
        with tempfile.TemporaryDirectory() as tmpdir:
            record = drv.StageRecord(tmpdir)
            with unittest.mock.patch.object(drv.subprocess, "Popen", return_value=fake), \
                 unittest.mock.patch.object(drv.os, "killpg", lambda pgid, sig: signals.append((pgid, sig))):
                with self.assertRaises(KeyboardInterrupt):
                    drv.run_child(["x"], os.path.join(tmpdir, "x.log"), 0, record, "sampling")
            # the whole process group of the stage gets SIGTERM; the record is removed
            self.assertEqual(signals, [(424242, signal.SIGTERM)])
            self.assertFalse(os.path.exists(record.path))

    def test_stage_record(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            record = drv.StageRecord(tmpdir)
            record.check()  # no record: nothing happens
            # run_child names the process group of the stage while it runs and clears it afterwards
            log = os.path.join(tmpdir, "x.log")
            probe = f"import os, sys; sys.exit(0 if os.path.exists({record.path!r}) else 3)"
            self.assertEqual(drv.run_child([sys.executable, "-c", probe], log, 0, record, "train"), 0)
            self.assertFalse(os.path.exists(record.path))
            self.assertEqual(drv.run_child([sys.executable, "-c", "import time; time.sleep(30)"], log, 1,
                                           record, "train"), -1)
            self.assertFalse(os.path.exists(record.path))
            # a live group running the recorded program (a stage that outlived a killed driver) is refused
            child = subprocess.Popen(["sleep", "60"], start_new_session=True)  # ps shows it as "sleep 60"
            try:
                record.write("sampling", child.pid, "/bin/sleep")
                with self.assertRaises(drv.DriverError) as cm:
                    record.check()
                self.assertIn("still running", str(cm.exception))
                self.assertIn(f"-{child.pid}", str(cm.exception))
                self.assertIn("after a reboot", str(cm.exception))
                # the same group id now runs something else (reused after a reboot): stale
                record.write("sampling", child.pid, "abics_sampling")
                with self.assertLogs("main", level="WARNING") as logs:
                    record.check()
                self.assertIn("not the sampling stage", logs.output[0])
                self.assertFalse(os.path.exists(record.path))
                record.write("sampling", child.pid, "sleep")
            finally:
                os.killpg(child.pid, 9)
                child.wait()
            # the group is gone: the record is stale and removed
            record.check()
            self.assertFalse(os.path.exists(record.path))
            # the driver's own group is never "a stage of an earlier driver"
            record.write("train", os.getpgrp(), "python")
            record.check()
            self.assertFalse(os.path.exists(record.path))
            # another host cannot be checked
            with open(record.path, "w") as f:
                f.write("some-other-host 12345 sampling 2026-01-01T00:00:00\n")
            with self.assertRaises(drv.DriverError) as cm:
                record.check()
            self.assertIn("some-other-host", str(cm.exception))
            # an unreadable record is removed
            with open(record.path, "w") as f:
                f.write("garbage\n")
            record.check()
            self.assertFalse(os.path.exists(record.path))

    def test_sigterm_becomes_an_exception(self):
        import signal

        previous = signal.getsignal(signal.SIGTERM)
        try:
            drv.install_signal_handlers()
            with self.assertRaises(drv.DriverInterrupted):
                os.kill(os.getpid(), signal.SIGTERM)
        finally:
            signal.signal(signal.SIGTERM, previous)

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
            # AL3 has no convergence.dat yet: the model trained on AL0..AL2 is not evaluated
            self.assertIsNone(meta["provenance"]["evaluated_in"])
            self.assertIn("has not been evaluated", meta["accuracy_note"])

            # the evaluation step (AL3 with statistics) describes the current model
            with open(os.path.join(tmpdir, "AL3", "convergence.dat"), "w") as f:
                f.write("# header\nall  2  0.05  0.05  0.05  0\n")
            drv.write_model_metadata(tmpdir, params_root, tomlfile, 3, driver_dir)
            with open(path) as f:
                meta = json.load(f)
            self.assertEqual(meta["provenance"]["evaluated_in"], "AL3")
            self.assertEqual(sorted(meta["accuracy"]), ["AL1", "AL2", "AL3"])
            self.assertAlmostEqual(meta["accuracy"]["AL3"]["rmse"], 0.05)
            self.assertIn("is evaluated in AL3", meta["accuracy_note"])
            self.assertEqual(meta["provenance"]["reference_configurations"], 12)  # training data unchanged


class TestRunDriverWithStubbedStages(unittest.TestCase):
    """The loop itself, with the three stages replaced by functions that write what the programs write."""

    params_root = {
        "sampling": {"nreplicas": 2, "solver": {"type": "aenet", "base_input_dir": "./baseinput"}},
        "mlref": {"nreplicas": 2, "ndata": 2, "solver": {"type": "user", "path": "./mock.py",
                                                         "run_scheme": "subprocess", "function": "m.f"}},
        "train": {"type": "aenet", "base_input_dir": "./aenet_train_input"},
        "config": {"supercell": [1, 1, 1], "base_structure": [], "defect_structure": []},
    }

    def stubs(self, tmpdir, calls, fail_sampling_at=None, fail_mlref_once=None, slow=None):
        def mlref(params_root):  # stands in for run_integrated_step
            n = drv.newest_mc_index(tmpdir)
            calls.append(("mlref", n))
            self.assertEqual(params_root["mlref"]["exec_mode"], "integrated")
            os.makedirs(os.path.join(tmpdir, f"AL{n}", "0"), exist_ok=True)
            if fail_mlref_once and n in fail_mlref_once:
                # interrupted: the directory exists, nothing is recorded
                fail_mlref_once.discard(n)
                raise drv.IntegratedStepError()
            if n > 0:  # the error statistics written from AL1 on
                with open(os.path.join(tmpdir, f"AL{n}", "convergence.dat"), "w") as f:
                    f.write(f"# header\nall  2  0.{n}  0.{n}  0.{n}  0\n")
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
            # freshly modified: refused, a sampling run may still be alive
            with self.assertRaises(drv.DriverError):
                self.run_loop(tmpdir, 1, calls)
            self.assertEqual(calls, [])
            make_old(os.path.join(tmpdir, "MC0"))
            self.assertEqual(self.run_loop(tmpdir, 1, calls), "done")
            self.assertEqual(calls, [("sampling", 0)])  # repeated, nothing else
            moved = [d for d in os.listdir(os.path.join(tmpdir, drv.DRIVER_DIR)) if d.startswith("MC0.interrupted-")]
            self.assertEqual(len(moved), 1)

    def test_stop_in_the_middle_of_an_iteration_describes_the_trained_model(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            calls = []
            # AL0 train MC0 AL1 train, then the sampling of MC1 fails
            with self.assertRaises(drv.DriverError):
                self.run_loop(tmpdir, 2, calls, fail_sampling_at=1)
            self.assertEqual(drv.read_progress(tmpdir), ["AL0", "train", "MC0", "AL1", "train"])
            calls.clear()
            self.assertEqual(self.run_loop(tmpdir, 1, calls), "done")  # target already met
            self.assertEqual(calls, [])
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "MC1")))  # the unrecorded MC1 is left alone
            with open(os.path.join(tmpdir, drv.DRIVER_DIR, drv.METADATA_FILE)) as f:
                meta = json.load(f)
            self.assertEqual(meta["provenance"]["iterations"], 1)
            self.assertEqual(meta["provenance"]["trained_on_steps"], 2)  # the model in baseinput used AL0 and AL1
            self.assertIn("AL0 .. AL1", meta["accuracy_note"])

    def test_max_iterations_below_one_is_refused(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(drv.DriverError):
                self.run_loop(tmpdir, 0, [])

    def run_with_al(self, tmpdir, al_section, max_iterations=2, evaluate=None, **stub_kwargs):
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
                return drv.run_driver(params_root, tomlfile, max_iterations, evaluate_final_model=evaluate), calls
        finally:
            os.chdir(cwd)

    def read_metadata(self, tmpdir):
        with open(os.path.join(tmpdir, drv.DRIVER_DIR, drv.METADATA_FILE)) as f:
            return json.load(f)

    def test_final_model_is_evaluated_by_one_more_reference_step(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state, calls = self.run_with_al(tmpdir, {"evaluate_final_model": True}, max_iterations=2)
            self.assertEqual(state, "done")
            self.assertEqual(calls, [("mlref", 0), ("train", 0), ("sampling", 0),
                                     ("mlref", 1), ("train", 1), ("sampling", 1), ("mlref", 2)])
            self.assertEqual(drv.read_progress(tmpdir), ["AL0", "train", "MC0", "AL1", "train", "MC1", "AL2"])
            meta = self.read_metadata(tmpdir)
            self.assertEqual(meta["provenance"]["iterations"], 2)
            self.assertEqual(meta["provenance"]["trained_on_steps"], 2)
            self.assertEqual(meta["provenance"]["evaluated_in"], "AL2")
            self.assertEqual(sorted(meta["accuracy"]), ["AL1", "AL2"])
            summary = [l for l in open(os.path.join(tmpdir, drv.SUMMARY_FILE)) if not l.startswith("#")]
            self.assertEqual([l.split()[0] for l in summary], ["1", "2"])
            # the evaluation exists: a rerun runs nothing
            state, calls = self.run_with_al(tmpdir, {"evaluate_final_model": True}, max_iterations=2)
            self.assertEqual((state, calls), ("done", []))
            # extending the loop trains on the evaluation data and evaluates the new model again
            state, calls = self.run_with_al(tmpdir, {"evaluate_final_model": True}, max_iterations=3)
            self.assertEqual(calls, [("train", 2), ("sampling", 2), ("mlref", 3)])
            self.assertEqual(self.read_metadata(tmpdir)["provenance"]["evaluated_in"], "AL3")

    def test_evaluation_can_be_requested_on_the_command_line(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state, calls = self.run_with_al(tmpdir, {}, max_iterations=1)
            self.assertEqual(calls, [("mlref", 0), ("train", 0), ("sampling", 0)])
            self.assertIsNone(self.read_metadata(tmpdir)["provenance"]["evaluated_in"])
            # the option given to run_driver (the --evaluate-final-model flag) adds the step afterwards
            state, calls = self.run_with_al(tmpdir, {}, max_iterations=1, evaluate=True)
            self.assertEqual((state, calls), ("done", [("mlref", 1)]))
            self.assertEqual(drv.read_progress(tmpdir), ["AL0", "train", "MC0", "AL1"])
            self.assertEqual(self.read_metadata(tmpdir)["provenance"]["evaluated_in"], "AL1")

    def test_interrupted_evaluation_is_resumed_and_a_lowered_target_skips_it(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            al = {"evaluate_final_model": True}
            # the evaluation AL1 is interrupted: AL1 exists, ALloop.progress still ends with MC0
            with self.assertRaises(drv.DriverError):
                self.run_with_al(tmpdir, al, max_iterations=1, fail_mlref_once={1})
            self.assertEqual(drv.read_progress(tmpdir), ["AL0", "train", "MC0"])
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "AL1")))
            # the rerun resumes the evaluation (run_integrated_step resumes an existing AL<n>)
            state, calls = self.run_with_al(tmpdir, al, max_iterations=1)
            self.assertEqual((state, calls), ("done", [("mlref", 1)]))
            self.assertEqual(drv.read_progress(tmpdir), ["AL0", "train", "MC0", "AL1"])
            # extended to 3 runs without evaluation, then a target below the recorded runs with the
            # option: the loop stops without spending reference calculations on an evaluation
            state, calls = self.run_with_al(tmpdir, {}, max_iterations=3)
            self.assertEqual(calls, [("train", 1), ("sampling", 1), ("mlref", 2), ("train", 2), ("sampling", 2)])
            self.assertEqual(drv.read_progress(tmpdir)[-1], "MC2")
            state, calls = self.run_with_al(tmpdir, al, max_iterations=2)
            self.assertEqual((state, calls), ("done", []))
            self.assertIsNone(self.read_metadata(tmpdir)["provenance"]["evaluated_in"])
            # with the target equal to the recorded runs the evaluation is added afterwards
            state, calls = self.run_with_al(tmpdir, al, max_iterations=3)
            self.assertEqual((state, calls), ("done", [("mlref", 3)]))

    def test_evaluation_respects_the_wall_time(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # the loop itself fits, the evaluation (an mlref stage) does not: stop with status "walltime"
            al = {"evaluate_final_model": True, "walltime": 100, "walltime_margin": 10,
                  "stage_estimates": {"mlref": 1, "train": 1, "sampling": 1}}
            state, calls = self.run_with_al(tmpdir, al, max_iterations=1)
            self.assertEqual((state, calls), ("done", [("mlref", 0), ("train", 0), ("sampling", 0), ("mlref", 1)]))
            with open(os.path.join(tmpdir, drv.DRIVER_DIR, drv.DURATIONS_FILE), "w") as f:
                json.dump({"mlref": 1000.0}, f)  # as if the reference-data step had taken 1000 s
            state, calls = self.run_with_al(tmpdir, al, max_iterations=2)
            self.assertEqual(state, "walltime")
            self.assertEqual(calls, [("train", 1), ("sampling", 1)])  # stopped before the evaluation AL2

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

    def test_driver_refuses_while_a_stage_of_a_killed_driver_runs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            driver_dir = os.path.join(tmpdir, drv.DRIVER_DIR)
            os.makedirs(driver_dir)
            child = subprocess.Popen(["sleep", "60"], start_new_session=True)
            try:
                drv.StageRecord(driver_dir).write("sampling", child.pid, "sleep")
                with self.assertRaises(drv.DriverError) as cm:
                    self.run_with_al(tmpdir, {})
                self.assertIn("still running", str(cm.exception))
                self.assertFalse(os.path.exists(os.path.join(driver_dir, drv.LOCK_FILE)))  # lock released
            finally:
                os.killpg(child.pid, 9)
                child.wait()
            state, calls = self.run_with_al(tmpdir, {}, max_iterations=1)  # stale record: removed, loop runs
            self.assertEqual(state, "done")
            self.assertEqual(calls, [("mlref", 0), ("train", 0), ("sampling", 0)])

    def test_lock_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            driver_dir = os.path.join(tmpdir, drv.DRIVER_DIR)
            os.makedirs(driver_dir)
            # a live process (this one, pretending to be abics_al) holds the lock
            with open(os.path.join(driver_dir, drv.LOCK_FILE), "w") as f:
                f.write(f"{os.uname().nodename} {os.getpid()} now\n")
            with unittest.mock.patch.object(drv, "_process_command", lambda pid: "python /venv/bin/abics_al input.toml"):
                with self.assertRaises(drv.DriverError):
                    drv.DriverLock(driver_dir).acquire()
            # ps cannot be read: refused as well (conservative)
            with unittest.mock.patch.object(drv, "_process_command", lambda pid: ""):
                with self.assertRaises(drv.DriverError):
                    drv.DriverLock(driver_dir).acquire()
            # the same pid reused by another program: the lock is stale
            with unittest.mock.patch.object(drv, "_process_command", lambda pid: "sleep 100"):
                lock = drv.DriverLock(driver_dir)
                with self.assertLogs("main", level="WARNING") as logs:
                    lock.acquire()
                self.assertIn("not abics_al", logs.output[0])
                lock.release()
            self.assertTrue(drv._process_command(os.getpid()))  # the real ps works here
            for command in ("/venv/bin/python /venv/bin/abics_al input.toml", "abics_al", "python -m abics.scripts.aldriver x"):
                self.assertTrue(drv._is_driver_command(command), command)
            for command in ("less al_driver/abics_al_1.log", "vim abics_al.lock", "sleep 100"):
                self.assertFalse(drv._is_driver_command(command), command)
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
