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

import random
import unittest

import numpy as np


from abics.mc import model, CanonicalMonteCarlo


class TestMC(unittest.TestCase):
    class model(model):
        def energy(self, config):
            return config[0] * config[1]

        def newconfig(self, config, dconfig):
            ret = config[:]
            ret[dconfig] *= -1.0
            return ret

        def trialstep(self, config, energy):
            dconfig = random.randint(0, 1)
            dE = self.energy(self.newconfig(config, dconfig)) - energy
            return dconfig, dE

    def test_mc(self):
        self.mc = CanonicalMonteCarlo(TestMC.model(), 1.0, [1.0,1.0])
        # self.mc = CanonicalMonteCarlo(TestMC.dmodel, 1.0, 1)
        N = 100
        S = 0.0
        S2 =0.0
        for i in range(N):
            E = self.mc.run(100, 1)[0]
            S += E
            S2 += E*E
        m = S/N
        er = np.sqrt((S2 - S*S/N)/(N-1)/N)
        self.assertTrue(np.isclose(-np.tanh(1.0), m, atol=3*er, rtol=0.0))
