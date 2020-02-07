import os
import unittest

import numpy as np

from pymatgen import Structure

from abics.applications.latgas_abinitio_interface.qe import QESolver


class TestQE(unittest.TestCase):
    def setUp(self):
        self.solver = QESolver('.')
        self.datadir = os.path.join(os.path.dirname(__file__), 'data', 'qe')

    # def test_get_results(self):
    #     res = self.solver.output.get_results(os.path.join(os.path.dirname(__file__), 'data', 'qe', 'ref'))
    #     ref = Structure.from_file(os.path.join(os.path.dirname(__file__), 'data', 'primitive.vasp'))
    #     res.structure.to('POSCAR', 'hoge.vasp')
    #     self.assertTrue(True)

    def test_input(self):
        self.solver.input.from_directory(os.path.join(self.datadir, 'baseinput'))
        A = np.ones((3,3)) - np.eye(3)
        r = np.array([0.0, 0.1, 0.2])
        st = Structure(A, ['H'], [r], coords_are_cartesian=False)
        self.solver.input.update_info_by_structure(st)
        self.assertTrue(np.allclose(self.solver.input.pwi.cell_parameters['cell'], A))
        self.assertTrue(np.allclose(self.solver.input.pwi.atomic_positions['positions'], r))
        self.solver.input.write_input(os.path.join(self.datadir, 'work'))



    def test_cl_algs(self):
        res = self.solver.input.cl_args(2, 4, 'work')
        self.assertEqual(res, ['-in', os.path.join('work', 'scf.in')])

