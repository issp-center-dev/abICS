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
import unittest

import numpy as np
from numpy.linalg import inv

from pymatgen import Structure

from abics.mc_mpi import RXParams
from abics.applications.latgas_abinitio_interface.model_setup import group
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)
from abics.applications.latgas_abinitio_interface.params import DFTParams


class TestInput(unittest.TestCase):
    def setUp(self):
        self.rootdir = os.path.dirname(__file__)
        self.datadir = os.path.join(self.rootdir, "data", "input")
        self.workdir = os.path.join(self.rootdir, "res", "input")

    def test_input(self):
        tomlfile = os.path.join(self.datadir, "input.toml")
        DFTParams.from_toml(tomlfile)

        rxparams = RXParams.from_toml(tomlfile)
        self.assertEqual(rxparams.nreplicas, 2)
        self.assertEqual(rxparams.nprocs_per_replica, 1)
        self.assertEqual(rxparams.kTstart, 1000.0)
        self.assertEqual(rxparams.kTend, 1200.0)
        self.assertEqual(rxparams.nsteps, 2)
        self.assertEqual(rxparams.RXtrial_frequency, 3)
        self.assertEqual(rxparams.sample_frequency, 4)
        self.assertEqual(rxparams.print_frequency, 5)
        self.assertEqual(rxparams.seed, 12345)

        configparams = DFTConfigParams.from_toml(tomlfile)
        spinel_config = defect_config(configparams)

        cellsize = np.array([2, 1, 1])
        self.assertTrue((spinel_config.cellsize == cellsize).all())

        A = np.zeros((3, 3))
        A[0, 0] = 2.0
        A[1, 1] = 1.0
        A[2, 2] = 0.5
        B = cellsize * A
        invB = inv(B)

        self.assertTrue(np.allclose(B, spinel_config.supercell))

        r = np.zeros((5, 3))
        r[1, :] = [0.5, 0.5, 0.5]
        r[2, :] = [0.0, 0.5, 0.5]
        r[3, :] = [0.5, 0.0, 0.5]
        r[4, :] = [0.5, 0.5, 0.0]
        st = Structure(A, ["O", "O", "Se", "Se", "Se"], r)

        seldyns = np.array(
            [
                [True, True, False],
                [True, True, False],
                [True, False, True],
                [True, False, True],
                [True, False, True],
                [True, False, True],
                [False, False, False],
                [False, False, False],
                [True, True, True],
                [True, True, True],
            ]
        )

        self.assertTrue(spinel_config.base_structure.matches(st))
        self.assertTrue(
            (spinel_config.base_structure.site_properties["seldyn"] == seldyns).all()
        )

        gs = [
            group("Al", ["Al"], coords=[[[0.0, 0.0, 0.0]]]),
            group("OH", ["O", "H"], coords=[[[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]]]),
        ]

        for i in range(2):
            g = spinel_config.defect_sublattices[0].groups[i]
            self.assertEqual(gs[i].name, g.name)
            self.assertEqual(gs[i].species, g.species)
            self.assertEqual(gs[i].natoms, g.natoms)
            self.assertTrue(np.allclose(np.dot(gs[i].coords, invB), g.coords))
