import importlib.util
import os
import sys
import tempfile
import unittest
import unittest.mock

import numpy as np
from pymatgen.core import Lattice, Structure

from abics.applications.latgas_abinitio_interface import activelearn as al
from abics.mlref import RefParams


def write_mc_replica(mc_replica_dir, energies, with_oxygen=True):
    """Write the per-replica output layout of DefaultObserver."""
    os.makedirs(mc_replica_dir)
    lattice = Lattice.cubic(4.0)
    with open(os.path.join(mc_replica_dir, "obs.dat"), "w") as f:
        f.write("# $1: step\n# $2: kT\n# $3: energy\n")
        for i, energy in enumerate(energies):
            f.write(f"{i}\t1000.0\t{energy}\n")
    np.save(
        os.path.join(mc_replica_dir, "obs_save.npy"),
        np.column_stack((np.array(energies), np.zeros(len(energies)))),
    )
    for i in range(len(energies)):
        # structure as seen by the sampling solver: O ignored
        st = Structure(lattice, ["Al", "Mg"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        st.to(fmt="POSCAR", filename=os.path.join(mc_replica_dir, f"structure.{i}.vasp"))
        species = ["Al", "Mg", "O"] if with_oxygen else ["Al", "Mg"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]][: len(species)]
        st_norel = Structure(lattice, species, coords)
        st_norel.to(
            fmt="POSCAR", filename=os.path.join(mc_replica_dir, f"structure_norel.{i}.vasp")
        )


class TestReferenceCandidates(unittest.TestCase):
    def test_default_candidates_follow_obs_dat(self):
        energies = [-1.0, -2.0, -3.0, -4.0]
        with tempfile.TemporaryDirectory() as tmpdir:
            mcdir = os.path.join(tmpdir, "MC0", "0")
            write_mc_replica(mcdir, energies)

            candidates = al.reference_candidates(mcdir, "RXMC")

            self.assertEqual([c.index for c in candidates], [0, 1, 2, 3])
            self.assertEqual([c.energy for c in candidates], energies)
            self.assertEqual([c.weight for c in candidates], [1.0] * 4)
            self.assertEqual(
                candidates[2].structure_file, os.path.join(mcdir, "structure.2.vasp")
            )
            self.assertEqual(
                candidates[2].structure_norel_file,
                os.path.join(mcdir, "structure_norel.2.vasp"),
            )

    def test_energies_come_from_the_line_of_the_structure(self):
        # print_frequency = 2 * sample_frequency: obs_save.npy has twice as
        # many rows as obs.dat, and structure.<i>.vasp belongs to line i of
        # obs.dat, not to row i of obs_save.npy.
        printed = [-1.0, -3.0, -5.0]
        saved = [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0]
        with tempfile.TemporaryDirectory() as tmpdir:
            mcdir = os.path.join(tmpdir, "MC0", "0")
            write_mc_replica(mcdir, printed)
            np.save(
                os.path.join(mcdir, "obs_save.npy"),
                np.column_stack((np.array(saved), np.zeros(len(saved)))),
            )

            candidates = al.reference_candidates(mcdir, "RXMC")

            self.assertEqual([c.index for c in candidates], [0, 1, 2])
            self.assertEqual([c.energy for c in candidates], printed)

    def test_energy_column_is_located_by_header_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, "obs.dat"), "w") as f:
                # a sampler with two parameters (e.g. kT and mu) before the observables
                f.write("# $1: step\n# $2: kT\n# $3: mu\n# $4: energy\n# $5: doi\n")
                f.write("0\t1000.0\t0.5\t-1.25\t0.1\n")
                f.write("1\t1000.0\t0.5\t-2.5\t0.2\n")
            self.assertEqual(al.read_obs_energies(tmpdir), [-1.25, -2.5])

            # without a header: step, temperature, energy
            with open(os.path.join(tmpdir, "obs.dat"), "w") as f:
                f.write("0\t1000.0\t-1.0\n1\t1000.0\t-2.0\n")
            self.assertEqual(al.read_obs_energies(tmpdir), [-1.0, -2.0])
            self.assertEqual(al.count_mc_samples(tmpdir), 2)

    def test_registered_reader_is_used_for_its_sampler_only(self):
        def reader(mc_replica_dir):
            return [al.ReferenceCandidate(7, "a.vasp", "b.vasp", 0.5, weight=2.0)]

        al.register_candidate_reader("testSampler", reader)
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                mcdir = os.path.join(tmpdir, "MC0", "0")
                write_mc_replica(mcdir, [-1.0, -2.0])

                custom = al.reference_candidates(mcdir, "testSampler")
                default = al.reference_candidates(mcdir, "PAMC")

            self.assertEqual([c.index for c in custom], [7])
            self.assertEqual(custom[0].weight, 2.0)
            self.assertEqual([c.index for c in default], [0, 1])
        finally:
            del al._candidate_readers["testSampler"]

    def test_refparams_reads_seed(self):
        params = RefParams.from_dict({"nreplicas": 2, "ndata": 5, "seed": 123})
        self.assertEqual(params.seed, 123)
        self.assertEqual(params.sampler, "linspace")
        self.assertEqual(RefParams.from_dict({"nreplicas": 2, "ndata": 5}).seed, 0)

    def test_select_candidates_linspace(self):
        candidates = [
            al.ReferenceCandidate(i, f"s{i}", f"n{i}", float(i)) for i in range(10)
        ]
        refparams = RefParams()
        refparams.ndata = 4
        refparams.sampler = "linspace"

        selected = al.select_candidates(candidates, refparams)

        self.assertEqual([c.index for c in selected], [0, 3, 6, 9])

    def test_select_candidates_rejects_ndata_larger_than_samples(self):
        candidates = [
            al.ReferenceCandidate(i, f"s{i}", f"n{i}", float(i)) for i in range(4)
        ]
        refparams = RefParams()
        refparams.ndata = 5
        for sampler in ("linspace", "random"):
            refparams.sampler = sampler
            with self.assertRaises(ValueError) as cm:
                al.select_candidates(candidates, refparams)
            self.assertIn("ndata = 5", str(cm.exception))
        refparams.ndata = 4
        refparams.sampler = "linspace"
        self.assertEqual(len(al.select_candidates(candidates, refparams)), 4)

    def test_select_candidates_random_is_a_sorted_subset(self):
        candidates = [
            al.ReferenceCandidate(i, f"s{i}", f"n{i}", float(i)) for i in range(10)
        ]
        refparams = RefParams()
        refparams.ndata = 4
        refparams.sampler = "random"

        selected = al.select_candidates(candidates, refparams)
        indices = [c.index for c in selected]

        self.assertEqual(len(set(indices)), 4)
        self.assertEqual(indices, sorted(indices))
        self.assertTrue(all(0 <= i < 10 for i in indices))

    def test_load_candidate_structures_restores_ignored_species(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mcdir = os.path.join(tmpdir, "MC0", "0")
            write_mc_replica(mcdir, [-1.0, -2.0])
            candidates = al.reference_candidates(mcdir, "RXMC")

            with_o = list(al.load_candidate_structures(candidates, ["O"]))
            without = list(al.load_candidate_structures(candidates, None))

        self.assertEqual(len(with_o), 2)
        self.assertEqual([site.species_string for site in with_o[0]], ["Al", "Mg", "O"])
        self.assertEqual([site.species_string for site in without[0]], ["Al", "Mg"])


class TestBookkeeping(unittest.TestCase):
    def test_stage_and_progress_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            aldir = os.path.join(tmpdir, "AL1")
            os.makedirs(aldir)
            self.assertEqual(al.read_stage(aldir), 0)
            al.record_stage(aldir, 0)
            self.assertEqual(al.read_stage(aldir), 1)
            al.record_stage(aldir, 1)
            self.assertEqual(al.read_stage(aldir), 2)

            al.record_progress(tmpdir, "AL0")
            al.record_progress(tmpdir, "train")
            self.assertEqual(al.last_progress(tmpdir), "train")

            al.write_rundirs(tmpdir, ["/a/input00/baseinput0", "/a/input01/baseinput0"])
            with open(os.path.join(tmpdir, "rundirs.txt")) as f:
                self.assertEqual(f.read(), "/a/input00/baseinput0\n/a/input01/baseinput0\n")

    def test_input_dir_format_and_indices(self):
        self.assertEqual(al.input_dir_format(10).format(3), "input03")
        self.assertEqual(al.input_dir_format(100).format(3), "input003")
        with tempfile.TemporaryDirectory() as tmpdir:
            for name in ("input00", "input08", "input17"):
                os.makedirs(os.path.join(tmpdir, name))
            self.assertEqual(al.list_input_indices(tmpdir), [0, 8, 17])

    def test_energy_corr_and_relax_max_format(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            corr = os.path.join(tmpdir, "energy_corr.dat")
            al.write_energy_corr(corr, [("NaN", -1.5, 0), (-2.0, -2.25, 8)])
            with open(corr) as f:
                self.assertEqual(f.read(), "NaN -1.5 0\n-2.0 -2.25 8\n")

            relax = os.path.join(tmpdir, "relax_max.dat")
            al.write_relax_max(relax, [(3, 0.05)])
            with open(relax) as f:
                self.assertEqual(f.read(), "3 \t 0.05\n")


class TestErrorStatistics(unittest.TestCase):
    def test_error_statistics_skips_nan_rows(self):
        corr = np.array(
            [
                [np.nan, -1.0, 0],
                [-2.0, -2.5, 8],  # error 0.5
                [-3.0, -2.0, 17],  # error 1.0
            ]
        )
        stats = al.error_statistics(corr)
        self.assertEqual(stats.n, 2)
        self.assertAlmostEqual(stats.rmse, np.sqrt((0.25 + 1.0) / 2))
        self.assertAlmostEqual(stats.mae, 0.75)
        self.assertAlmostEqual(stats.max_error, 1.0)
        self.assertEqual(stats.max_error_sample, 17)

    def test_error_statistics_without_valid_rows(self):
        stats = al.error_statistics(np.array([[np.nan, -1.0, 0]]))
        self.assertEqual(stats.n, 0)
        self.assertTrue(np.isnan(stats.rmse))
        self.assertEqual(stats.max_error_sample, -1)
        self.assertEqual(al.error_statistics(np.zeros((0, 3))).n, 0)

    def test_al_step_statistics_and_convergence_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            aldir = os.path.join(tmpdir, "AL1")
            os.makedirs(os.path.join(aldir, "0"))
            os.makedirs(os.path.join(aldir, "1"))
            al.write_energy_corr(
                os.path.join(aldir, "0", "energy_corr.dat"), [(-1.0, -1.5, 0), (-2.0, -2.0, 8)]
            )
            al.write_energy_corr(
                os.path.join(aldir, "1", "energy_corr.dat"), [(-1.0, -1.0, 0), (-3.0, -2.0, 8)]
            )

            rows = al.read_al_step_rows(aldir, 2)
            self.assertEqual([r.shape for r in rows], [(2, 3), (2, 3)])
            per_replica, overall = al.al_step_error_statistics(rows)
            # the same from in-memory rows as gathered by abics_mlref (lists of tuples)
            per_replica_mem, overall_mem = al.al_step_error_statistics(
                [[(-1.0, -1.5, 0), (-2.0, -2.0, 8)], [(-1.0, -1.0, 0), (-3.0, -2.0, 8)]]
            )
            self.assertEqual(per_replica_mem, per_replica)
            self.assertEqual(overall_mem, overall)
            self.assertEqual([r for r, _ in per_replica], [0, 1])
            self.assertAlmostEqual(per_replica[0][1].max_error, 0.5)
            self.assertEqual(per_replica[0][1].max_error_sample, 0)
            self.assertAlmostEqual(per_replica[1][1].max_error, 1.0)
            self.assertEqual(per_replica[1][1].max_error_sample, 8)
            self.assertEqual(overall.n, 4)
            self.assertAlmostEqual(overall.mae, (0.5 + 0.0 + 0.0 + 1.0) / 4)
            self.assertAlmostEqual(overall.rmse, np.sqrt((0.25 + 1.0) / 4))

            lines = al.format_convergence("AL1", per_replica, overall)
            path = os.path.join(aldir, "convergence.dat")
            al.write_convergence(path, lines)
            with open(path) as f:
                content = f.read().splitlines()
            data = [line.split() for line in content if not line.startswith("#")]
            self.assertEqual([row[0] for row in data], ["0", "1", "all"])
            self.assertEqual(data[2][1], "4")
            self.assertAlmostEqual(float(data[2][3]), 0.375)


class TestParityPlot(unittest.TestCase):
    rows = [
        [(-1.0, -1.5, 0), (-2.0, -2.0, 8)],
        [(float("nan"), -1.0, 0), (-3.0, -2.0, 8)],
    ]

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib not installed")
    def test_parity_plot_is_written(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "parity.png")
            self.assertTrue(al.write_parity_plot(path, "AL1", self.rows))
            self.assertGreater(os.path.getsize(path), 0)
            with open(path, "rb") as f:
                self.assertEqual(f.read(8), b"\x89PNG\r\n\x1a\n")

    def test_parity_plot_is_skipped_without_matplotlib(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "parity.png")
            with unittest.mock.patch.dict(sys.modules, {"matplotlib": None}):
                self.assertFalse(al.write_parity_plot(path, "AL1", self.rows))
            self.assertFalse(os.path.exists(path))


if __name__ == "__main__":
    unittest.main()
