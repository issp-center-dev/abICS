import tempfile
import unittest
from pathlib import Path
from typing import Dict
from unittest.mock import patch

try:
    from abics.applications.latgas_abinitio_interface.nequip import NequipSolver
except ModuleNotFoundError as exc:
    if exc.name not in ("torch", "nequip"):
        raise
    NequipSolver = None


if NequipSolver is not None:
    import torch
    import yaml
    from pymatgen.core import Lattice, Structure

    class TypeEnergyModel(torch.nn.Module):
        def forward(self, data: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
            return {"total_energy": data["atom_types"].to(torch.float64).sum()}


@unittest.skipIf(
    NequipSolver is None, "Optional torch/nequip dependencies are not installed"
)
class TestNequipSolverInput(unittest.TestCase):

    def test_deployed_type_names_take_precedence_over_chemical_symbols(self):
        """
        The atom type ordering used by a deployed NequIP model must come from
        the deployed model metadata, not from input.yaml chemical_symbols.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            base_input_dir = Path(tmpdir)

            # input.yaml intentionally has a different ordering.
            (base_input_dir / "input.yaml").write_text(
                "chemical_symbols:\n"
                "  - Cu\n"
                "  - V\n"
                "  - Sn\n"
                "r_max: 8.0\n"
            )

            # torch.jit.load only needs the file path to exist for this mock.
            (base_input_dir / "deployed.pth").touch()

            dummy_model = object()

            def fake_jit_load(path, _extra_files=None):
                self.assertEqual(Path(path), base_input_dir / "deployed.pth")

                if _extra_files is not None:
                    _extra_files["type_names"] = b"V Cu Sn"
                    _extra_files["r_max"] = b"8.0"

                return dummy_model

            with patch(
                "abics.applications.latgas_abinitio_interface.nequip.torch.jit.load",
                side_effect=fake_jit_load,
            ):
                solver_input = NequipSolver.Input(ignore_species=None)
                solver_input.from_directory(str(base_input_dir))

            self.assertIs(solver_input.model, dummy_model)
            self.assertEqual(solver_input.element_list, ["V", "Cu", "Sn"])
            self.assertEqual(solver_input.r_max, 8.0)

    def test_custom_type_names_use_deployed_chemical_mapping(self):
        # Both mapping directions are supported, including YAML string keys.
        for mapping in (
            {"chemical_symbol_to_type": {"H": 1, "O": 0}},
            {"type_to_chemical_symbol": {"0": "O", "1": "H"}},
        ):
            with self.subTest(mapping=mapping), tempfile.TemporaryDirectory() as tmpdir:
                model_path = Path(tmpdir) / "deployed.pth"
                model = torch.jit.script(TypeEnergyModel())
                config = dict(mapping, _jit_fusion_strategy=[("DYNAMIC", 3)])
                torch.jit.save(
                    model, str(model_path), _extra_files={
                        "r_max": "8.0",
                        "type_names": "oxygen hydrogen",
                        "config": yaml.dump(config),
                    },
                )
                # A local input must not override the deployed configuration.
                (Path(tmpdir) / "input.yaml").write_text("chemical_symbols: [H, O]\n")
                solver = NequipSolver(ignore_species=None)
                solver.input.from_directory(tmpdir)
                self.assertEqual(solver.input.element_list, ["O", "H"])
                structure = Structure(
                    Lattice.cubic(10),
                    ["H", "H", "O"],
                    [[0, 0, 0], [0.1, 0, 0], [0, 0.1, 0]],
                )
                solver.input.update_info_by_structure(structure)
                solver.calculate_energy(None, tmpdir)
                # Two hydrogens of type 1 and one oxygen of type 0.
                self.assertEqual(solver.output.ene, 2.0)

    def test_custom_type_names_with_chemical_symbols_only(self):
        # Legacy configs may only specify chemical_symbols. NequIP sorts them
        # by atomic number, regardless of their order in the YAML file.
        with tempfile.TemporaryDirectory() as tmpdir:
            torch.jit.save(
                torch.jit.script(TypeEnergyModel()),
                str(Path(tmpdir) / "deployed.pth"),
                _extra_files={
                    "r_max": "8.0",
                    "type_names": "hydrogen oxygen",
                    "config": "chemical_symbols: [O, H]\n",
                },
            )
            solver_input = NequipSolver.Input()
            solver_input.from_directory(tmpdir)
            self.assertEqual(solver_input.element_list, ["H", "O"])

    def test_string_metadata_without_local_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            def fake_jit_load(path, _extra_files):
                _extra_files.update(type_names="V Cu Sn", r_max="8.0", config="")
                return object()

            with patch(
                "abics.applications.latgas_abinitio_interface.nequip.torch.jit.load",
                side_effect=fake_jit_load,
            ):
                solver_input = NequipSolver.Input()
                solver_input.from_directory(tmpdir)
            self.assertEqual(solver_input.element_list, ["V", "Cu", "Sn"])
            self.assertEqual(solver_input.r_max, 8.0)

    def test_deployed_allegro_config_ignores_unrelated_python_objects(self):
        from e3nn.o3 import Irreps

        # Allegro writes Irreps objects into its deployed YAML configuration.
        # These are unrelated to the chemical mapping and must not be loaded.
        for mapping, type_names, elements in (
            ({"chemical_symbol_to_type": {"H": 1, "O": 0}},
             "oxygen hydrogen", ["O", "H"]),
            ({"type_to_chemical_symbol": {0: "O", 1: "H"}},
             "oxygen hydrogen", ["O", "H"]),
            ({"chemical_symbols": ["O", "H"]},
             "hydrogen oxygen", ["H", "O"]),
            ({}, "H O", ["H", "O"]),
        ):
            with self.subTest(mapping=mapping), tempfile.TemporaryDirectory() as tmpdir:
                config = dict(
                    mapping,
                    tensor_track_allowed_irreps=Irreps("0e + 1o"),
                    _jit_fusion_strategy=[("DYNAMIC", 3)],
                )
                config_yaml = yaml.dump(config)
                self.assertIn("!!python/object/new:e3nn.o3._irreps.Irreps", config_yaml)
                torch.jit.save(
                    torch.jit.script(TypeEnergyModel()),
                    str(Path(tmpdir) / "deployed.pth"),
                    _extra_files={
                        "r_max": "8.0",
                        "type_names": type_names,
                        "config": config_yaml,
                    },
                )
                solver_input = NequipSolver.Input()
                solver_input.from_directory(tmpdir)
                self.assertEqual(solver_input.element_list, elements)
                self.assertEqual(solver_input.r_max, 8.0)


if __name__ == "__main__":
    unittest.main()
