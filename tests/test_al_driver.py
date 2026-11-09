import json
import os
import sys
import tempfile
import unittest

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

    def test_inconsistent_state_is_refused(self):
        with self.assertRaises(drv.DriverError):
            drv.next_action([], 1, 2)  # MC0 exists without a progress file
        with self.assertRaises(drv.DriverError):
            drv.next_action(["AL0"], 1, 2)  # MC0 exists but the progress says AL0 only
        with self.assertRaises(drv.DriverError):
            drv.next_action(["AL0", "train", "MC0"], 2, 3)  # MC1 exists but is not recorded
        with self.assertRaises(drv.DriverError):
            drv.next_action(["AL0", "train"], 1, 2)  # 'train' after AL0 but MC0 already exists
        with self.assertRaises(drv.DriverError):
            drv.next_action(["something"], 0, 2)


class TestParamsAndClock(unittest.TestCase):
    def test_driver_params(self):
        p = drv.DriverParams.from_dict({"max_iterations": 4, "walltime": 3600, "resubmit_command": "sbatch job.sh"})
        self.assertEqual((p.max_iterations, p.walltime, p.resubmit_command), (4, 3600.0, "sbatch job.sh"))
        self.assertIn("{input}", p.train_command)
        self.assertIn("{nreplicas}", p.sampling_command)
        with self.assertRaises(ValueError):
            drv.DriverParams.from_dict({"max_iterations": 0})
        with self.assertRaises(ValueError):
            drv.DriverParams.from_dict({"walltime": -1})

    def test_wall_clock(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            params = drv.DriverParams.from_dict({"walltime": 100, "walltime_margin": 10})
            clock = drv.WallClock(params, tmpdir)
            ok, _ = clock.can_start("mlref")  # nothing known yet: only the margin counts
            self.assertTrue(ok)
            clock.record("mlref", 95.0)  # the last mlref took 95 s: it would not fit any more
            ok, reason = clock.can_start("mlref")
            self.assertFalse(ok)
            self.assertIn("95 s", reason)
            ok, _ = clock.can_start("train")  # unknown duration: allowed
            self.assertTrue(ok)
            # durations persist across runs
            with open(os.path.join(tmpdir, drv.DURATIONS_FILE)) as f:
                self.assertEqual(json.load(f), {"mlref": 95.0})
            self.assertEqual(drv.WallClock(params, tmpdir).durations, {"mlref": 95.0})

            unlimited = drv.WallClock(drv.DriverParams.from_dict({}), tmpdir)
            self.assertTrue(unlimited.can_start("mlref")[0])


class TestStages(unittest.TestCase):
    def test_run_child(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log = os.path.join(tmpdir, "x.log")
            self.assertEqual(drv.run_child(f"{sys.executable} -c 'print(42)'", log, 0), 0)
            with open(log) as f:
                self.assertEqual(f.read().strip(), "42")
            self.assertEqual(drv.run_child(f"{sys.executable} -c 'import sys; sys.exit(3)'", log, 0), 3)
            self.assertEqual(drv.run_child(f"{sys.executable} -c 'import time; time.sleep(30)'", log, 1), -1)
            self.assertEqual(drv.run_child("no-such-command-abics x", log, 0), 127)

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
            base = os.path.join(tmpdir, "aenet_train_input")
            params_root["train"]["base_input_dir"] = base
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

            lines = drv.write_summary(tmpdir)
            data = [l.split() for l in lines if not l.startswith("#")]
            self.assertEqual(data, [["1", "4", "0.158", "0.15", "0.2"], ["2", "4", "0.079", "0.075", "0.1"]])
            self.assertTrue(os.path.isfile(os.path.join(tmpdir, drv.SUMMARY_FILE)))

            self.assertEqual(drv.count_training_data(tmpdir), 12)

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
            self.assertEqual(meta["provenance"]["training_configurations"], 12)
            self.assertEqual(meta["provenance"]["reference_solver"]["type"], "qe")
            self.assertEqual(meta["provenance"]["mlref_seed"], 7)
            self.assertIsNone(meta["provenance"]["base_model"])
            self.assertEqual(meta["scope"]["species"], ["Al", "Mg", "O"])
            self.assertEqual(meta["scope"]["kT_range"], [1000.0, 1500.0])
            self.assertEqual(meta["model_files"], "./baseinput")
            self.assertEqual(sorted(meta["accuracy"]), ["AL1", "AL2"])
            self.assertAlmostEqual(meta["accuracy"]["AL2"]["rmse"], 0.079)
            self.assertEqual(len(meta["input_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
