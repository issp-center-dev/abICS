# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import os
import shutil
import unittest

import numpy as np

from pymatgen import Structure
from pymatgen.io.vasp.inputs import VaspInput

from abics.applications.latgas_abinitio_interface.vasp import VASPSolver


class TestVASP(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rootdir = os.path.dirname(__file__)
        workdir = os.path.join(rootdir, "res", "vasp")
        if os.path.exists(workdir):
            shutil.rmtree(workdir)
        os.makedirs(workdir)

    def setUp(self):
        self.solver = VASPSolver(".")
        self.rootdir = os.path.dirname(__file__)
        self.datadir = os.path.join(self.rootdir, "data", "vasp")
        self.workdir = os.path.join(self.rootdir, "res", "vasp")

    def test_get_results(self):
        res = self.solver.output.get_results(os.path.join(self.datadir, "output"))
        res.structure.to("POSCAR", os.path.join(self.workdir, "pos.vasp"))
        ref = Structure.from_file(os.path.join(self.datadir, "..", "pos.vasp"))
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
            site_properties={"seldyn": [[True, True, False], [True, False, True]]},
        )
        self.solver.input.update_info_by_structure(st)
        self.solver.input.write_input(self.workdir)

        res = VaspInput.from_directory(self.workdir)
        ref = VaspInput.from_directory(os.path.join(self.datadir, "input"))

        self.assertEqual(res["INCAR"], ref["INCAR"])
        self.assertTrue(res["POSCAR"].structure.matches(ref["POSCAR"].structure))
        self.assertEqual(
            res["POSCAR"].structure.site_properties,
            ref["POSCAR"].structure.site_properties,
        )

    def test_cl_algs(self):
        nprocs = 2
        nthreads = 4
        workdir = "work"
        res = self.solver.input.cl_args(nprocs, nthreads, workdir)
        self.assertEqual(res, [workdir])

