"""Check the example observer and reference-output failure handling.

TorchSim tests are optional: normal abICS installations need no GPU packages.
"""
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from pymatgen.core import Structure


EXAMPLE = Path(__file__).resolve().parents[1] / 'examples' / 'active_learning_umlp'


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, EXAMPLE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestUmlpObserver(unittest.TestCase):
    def test_doi_and_import_from_another_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            previous = os.getcwd()
            try:
                os.chdir(temporary)
                observer = load_module('umlp_observer', 'calc_DOI.py')
            finally:
                os.chdir(previous)
        structure = Structure.from_file(EXAMPLE / 'MgAl2O4.vasp')
        self.assertEqual(observer.calc_DOI(structure), 0.0)
        mg = structure.indices_from_symbol('Mg')
        al = structure.indices_from_symbol('Al')
        for i, j in zip(mg, al):
            structure.replace(i, 'Al')
            structure.replace(j, 'Mg')
        self.assertEqual(observer.calc_DOI(structure), 1.0)


HAS_TORCHSIM = (
    sys.version_info >= (3, 12)
    and importlib.util.find_spec('torch_sim') is not None
    and importlib.util.find_spec('orb_models') is not None
)


@unittest.skipUnless(HAS_TORCHSIM, 'requires optional TorchSim example dependencies')
class TestUmlpReferenceOutputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner = load_module('umlp_torchsim_runner', 'run_umlp_torchsim.py')

    def setUp(self):
        import ase.io
        import torch
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.dirs = [self.root / 'first', self.root / 'second']
        self.atoms = []
        for i, directory in enumerate(self.dirs):
            directory.mkdir()
            atoms = ase.io.read(EXAMPLE / 'MgAl2O4.vasp')
            atoms.positions += i * 0.03
            ase.io.write(directory / 'structure.vasp', atoms, format='vasp')
            self.atoms.append(atoms)
        self.initial = [(p / 'structure.vasp').read_bytes() for p in self.dirs]
        self.rundirs = self.root / 'rundirs.txt'
        self.rundirs.write_text(''.join(str(p) + '\n' for p in self.dirs))
        self.summary = self.root / 'summary.json'
        self.argv = ['run_umlp_torchsim.py', '--rundirs', str(self.rundirs),
                     '--summary', str(self.summary), '--device', 'cpu']
        self.model = SimpleNamespace(device=torch.device('cpu'), max_batch_seen=2,
                                     forward_calls=4)

    def run_with_results(self, converged=(True, True), energy=-10.0):
        final = [a.copy() for a in self.atoms]
        for a in final:
            a.positions += 0.01
        records = [{'converged': ok, 'energy_eV': energy - i,
                    'atomic_fmax': 0.01, 'cell_fmax': 0.02}
                   for i, ok in enumerate(converged)]
        with patch.object(sys, 'argv', self.argv), \
             patch.object(self.runner, 'load_model', return_value=self.model), \
             patch.object(self.runner, 'relax_batch', return_value=(final, records)), \
             contextlib.redirect_stdout(io.StringIO()):
            self.runner.main()

    def test_converged_results_keep_input_order(self):
        import ase.io
        import numpy as np
        self.run_with_results()
        for i, directory in enumerate(self.dirs):
            self.assertEqual(float((directory / 'energy.dat').read_text()), -10.0 - i)
            actual = ase.io.read(directory / 'structure.vasp')
            np.testing.assert_allclose(actual.positions, self.atoms[i].positions + .01)
            original = ase.io.read(directory / 'structure_norel.vasp')
            np.testing.assert_allclose(original.positions, self.atoms[i].positions)
        self.assertEqual(len(json.loads(self.summary.read_text())['structures']), 2)

    def test_unconverged_structure_gets_no_reference_energy(self):
        with self.assertRaisesRegex(SystemExit, '1 unconverged references'):
            self.run_with_results(converged=(True, False))
        self.assertTrue((self.dirs[0] / 'energy.dat').exists())
        self.assertFalse((self.dirs[1] / 'energy.dat').exists())
        self.assertEqual((self.dirs[1] / 'structure.vasp').read_bytes(), self.initial[1])
        report = json.loads((self.dirs[1] / 'torchsim_relaxation.json').read_text())
        self.assertFalse(report['converged'])
        self.assertTrue((self.dirs[1] / 'structure_unconverged.vasp').exists())

    def test_nonfinite_energy_cannot_be_collected(self):
        with self.assertRaisesRegex(RuntimeError, 'Non-finite'):
            self.run_with_results(energy=float('nan'))
        for directory, original in zip(self.dirs, self.initial):
            self.assertFalse((directory / 'energy.dat').exists())
            self.assertEqual((directory / 'structure.vasp').read_bytes(), original)

    def test_existing_outputs_are_not_overwritten(self):
        energy = self.dirs[1] / 'energy.dat'
        energy.write_text('123.0\n')
        with patch.object(sys, 'argv', self.argv), \
             patch.object(self.runner, 'load_model') as loader, \
             contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            self.runner.main()
        loader.assert_not_called()
        self.assertEqual(energy.read_text(), '123.0\n')
        self.assertFalse((self.dirs[0] / 'energy.dat').exists())

    def test_nonfinite_tolerance_is_rejected_before_loading(self):
        with patch.object(sys, 'argv', self.argv + ['--fmax', 'nan']), \
             patch.object(self.runner, 'load_model') as loader, \
             contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            self.runner.main()
        loader.assert_not_called()
