import unittest
from unittest.mock import patch

from abics.applications.latgas_abinitio_interface.nequip import NequipSolver


class TestNequipSolverInput(unittest.TestCase):

    def test_deployed_type_names_take_precedence_over_chemical_symbols(self):
        """
        The atom type ordering used by a deployed NequIP model must come from
        the deployed model metadata, not from input.yaml chemical_symbols.
        """
        import tempfile
        from pathlib import Path

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


if __name__ == "__main__":
    unittest.main()
