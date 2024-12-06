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

from pymatgen.core import Structure

from abics.applications.latgas_abinitio_interface.base_solver import create_solver
from abics.applications.latgas_abinitio_interface.params import DFTParams
from abics.applications.latgas_abinitio_interface.aenet_pylammps import (
    AenetPyLammpsSolver,
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
