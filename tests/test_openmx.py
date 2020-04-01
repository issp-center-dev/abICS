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

from abics.applications.latgas_abinitio_interface.openmx import OpenMXInputFile, OpenMXSolver


class TestOpenMX(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rootdir = os.path.dirname(__file__)
        workdir = os.path.join(rootdir, "res", "openmx")
        if os.path.exists(workdir):
            shutil.rmtree(workdir)
        os.makedirs(workdir)

    def setUp(self):
        self.rootdir = os.path.dirname(__file__)
        self.datadir = os.path.join(self.rootdir, "data", "openmx")
        self.workdir = os.path.join(self.rootdir, "res", "openmx")
        self.solver = OpenMXSolver(os.path.join(self.datadir, "bin", "openmx.dummy"))

    def test_get_results(self):
        self.solver.input.from_directory(os.path.join(self.datadir, "baseinput"))
        res = self.solver.output.get_results(os.path.join(self.datadir, "output"))
        res.structure.to("POSCAR", os.path.join(self.workdir, "pos.vasp"))
        ref = Structure.from_file(os.path.join(self.datadir, "..", "pos.vasp"))
        ref_energy = -119.28626359359154
        self.assertTrue(res.structure.matches(ref))
        self.assertTrue(np.isclose(res.energy, ref_energy))

    def test_input(self):
        self.solver.input.from_directory(os.path.join(self.datadir, "baseinput"))
        A = 4.0*np.eye(3)
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
        res = OpenMXInputFile(os.path.join(self.workdir, "Al.dat"))
        ref = OpenMXInputFile(os.path.join(self.datadir, "input", "Al.dat"))
        ref["System.CurrrentDirectory"] = [self.workdir + "/"]

        self.assertEqual(res, ref)


    def test_cl_algs(self):
        self.solver.input.from_directory(os.path.join(self.datadir, "baseinput"))
        nprocs = 2
        nthreads = 4
        workdir = "work"
        res = self.solver.input.cl_args(nprocs, nthreads, workdir)
        self.assertEqual(res, [os.path.join(workdir, "Al.dat"), "-nt", str(nthreads)])

