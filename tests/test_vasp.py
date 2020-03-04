import os
import unittest

import numpy as np

from pymatgen import Structure
from pymatgen.io.vasp.inputs import VaspInput

from abics.applications.latgas_abinitio_interface.vasp import VASPSolver


class TestVASP(unittest.TestCase):
    def setUp(self):
        self.solver = VASPSolver(".")
        self.rootdir = os.path.dirname(__file__)
        self.datadir = os.path.join(self.rootdir, "data", "vasp")
        self.workdir = os.path.join(self.rootdir, "res", "vasp")

    def test_get_results(self):
        res = self.solver.output.get_results(os.path.join(self.datadir, "output"))
        res.structure.to("POSCAR", os.path.join(self.workdir, "pos.vasp"))
        ref = Structure.from_file(os.path.join(self.datadir, "pos.vasp"))
        ref_energy = 0.54083824
        self.assertTrue(np.isclose(res.energy, ref_energy))
        self.assertTrue(res.structure.matches(ref))

    def test_input(self):
        self.solver.input.from_directory(os.path.join(self.datadir, "baseinput"))
        A = 4.0 * np.eye(3)
        r = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        st = Structure(
            A,
            ["Al", "Al"],
            r,
            coords_are_cartesian=False,
            site_properties={"seldyn": [[True, True, True], [True, True, True]]},
        )
        self.solver.input.update_info_by_structure(st)
        self.solver.input.write_input(self.workdir)

        res = VaspInput.from_directory(self.workdir)
        ref = VaspInput.from_directory(os.path.join(self.datadir, "input"))

        self.assertEqual(res["INCAR"], ref["INCAR"])
        self.assertTrue(res["POSCAR"].structure.matches(ref["POSCAR"].structure))

    def test_cl_algs(self):
        nprocs = 2
        nthreads = 4
        workdir = "work"
        res = self.solver.input.cl_args(nprocs, nthreads, workdir)
        self.assertEqual(res, [workdir])

