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
import tempfile
import unittest

import numpy as np

from pymatgen.core import Structure

from abics.applications.latgas_abinitio_interface.base_solver import create_solver
from abics.applications.latgas_abinitio_interface.params import DFTParams
from abics.applications.latgas_abinitio_interface.aenet_pylammps import (
    AenetPyLammpsSolver,
    get_lammps_species_map,
)


class TestAenetPyLammps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rootdir = os.path.dirname(__file__)
        workdir = os.path.join(rootdir, "res", "aenet")
        if os.path.exists(workdir):
            shutil.rmtree(workdir)
        os.makedirs(workdir)

    def setUp(self):
        params = DFTParams.from_dict({"type": "aenetPyLammps", "path": "."})

        self.imported = False
        try:
            self.solver = create_solver(params.solver, params)
            self.imported = True
        except ImportError:
            self.imported = False

        self.rootdir = os.path.dirname(__file__)
        self.datadir = os.path.join(self.rootdir, "data", "aenet")
        self.workdir = os.path.join(self.rootdir, "res", "aenet")

    def test_create_solver(self):
        if self.imported:
            self.assertIsInstance(self.solver, AenetPyLammpsSolver)

    def test_input_reads_species_order_from_in_lammps(self):
        input_mgr = AenetPyLammpsSolver.Input()
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, "in.lammps"), "w") as f:
                f.write("pair_style aenet\n")
                f.write("pair_coeff * * v00 Al Mg 15t-15t.nn Al Mg\n")
            input_mgr.from_directory(tmpdir)

        self.assertEqual(input_mgr.lammps_species, ["Al", "Mg"])

    def test_species_map_uses_in_lammps_order(self):
        st = Structure(
            np.eye(3),
            ["Mg"],
            [[0.0, 0.0, 0.0]],
            coords_are_cartesian=False,
        )
        spec_dict, nspec = get_lammps_species_map(st, ["Al", "Mg"])

        self.assertEqual(spec_dict["Mg"], 2)
        self.assertEqual(nspec, 2)
