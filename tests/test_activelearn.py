import os
import tempfile
import unittest

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
    def test_default_candidates_follow_obs_dat_and_obs_save(self):
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

    def test_select_candidates_linspace(self):
        candidates = [
            al.ReferenceCandidate(i, f"s{i}", f"n{i}", float(i)) for i in range(10)
        ]
        refparams = RefParams()
        refparams.ndata = 4
        refparams.sampler = "linspace"

        selected = al.select_candidates(candidates, refparams)

        self.assertEqual([c.index for c in selected], [0, 3, 6, 9])

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


if __name__ == "__main__":
    unittest.main()
