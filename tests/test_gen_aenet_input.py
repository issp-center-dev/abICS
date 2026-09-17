import os
import tempfile
import unittest

from abics.scripts.gen_aenet_input import get_species, main_impl


class TestGenAenetInput(unittest.TestCase):
    def test_get_species_includes_all_configured_groups(self):
        params_root = {
            "config": {
                "unitcell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "supercell": [1, 1, 1],
                "base_structure": [{}],
                "defect_structure": [
                    {
                        "coords": [[0, 0, 0]],
                        "groups": [
                            {"name": "Al", "num": 1},
                            {"name": "Mg", "num": 0},
                        ],
                    }
                ],
                "chemical_potential": [{"species": "Mg", "mu": 0.1}],
            },
            "train": {"ignore_species": []},
        }

        self.assertEqual(get_species(params_root), ["Al", "Mg"])

    def test_main_impl_uses_sampling_ignore_species_for_lammps(self):
        params_root = {
            "sampling": {
                "solver": {"type": "aenetPyLammps", "ignore_species": ["O"]},
            },
            "train": {"type": "aenet", "ignore_species": [], "base_input_dir": "unused"},
            "config": {
                "unitcell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "supercell": [1, 1, 1],
                "base_structure": [{"type": "O", "coords": [[0, 0, 0]]}],
                "defect_structure": [
                    {
                        "coords": [[0.5, 0.5, 0.5]],
                        "groups": [
                            {"name": "Al", "num": 1},
                            {"name": "Mg", "num": 0},
                        ],
                    }
                ],
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            main_impl(params_root, tmpdir, force=False)
            with open(os.path.join(tmpdir, "predict", "in.lammps")) as f:
                self.assertIn("pair_coeff      * * v00 Al Mg 15t-15t.nn Al Mg", f.read())

    def test_main_impl_rejects_untrained_lammps_species(self):
        params_root = {
            "sampling": {
                "solver": {"type": "aenetPyLammps", "ignore_species": []},
            },
            "train": {"type": "aenet", "ignore_species": ["Mg"], "base_input_dir": "unused"},
            "config": {
                "unitcell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "supercell": [1, 1, 1],
                "base_structure": [{}],
                "defect_structure": [
                    {
                        "coords": [[0.5, 0.5, 0.5]],
                        "groups": [
                            {"name": "Al", "num": 1},
                            {"name": "Mg", "num": 0},
                        ],
                    }
                ],
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "untrained species"):
                main_impl(params_root, tmpdir, force=False)
