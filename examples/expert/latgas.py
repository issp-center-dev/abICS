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

import numpy as np
import random as rand
import sys

from abics.mc import model, CanonicalMonteCarlo, binning, observer_base


class latticegas(model):
    """This class defines the 2D lattice gas model"""

    model_name = "latticegas"

    def __init__(self, Eads, J=0, mu=0):
        self.J = J
        self.Eads = Eads
        self.mu = mu

    def energy(self, latgas_config):
        """ Calculate total energy of the 2D lattice gas model"""
        e = 0.0
        config = latgas_config.config  # config should be 2D numpy array
        for i in range(latgas_config.lenX):
            for j in range(latgas_config.lenY):
                e += self.J * config[i, j] * (config[i - 1, j] + config[i, j - 1])
        e += config.sum() * (self.Eads - self.mu)
        return e

    def density(self, latgas_config):
        """Calculate number of particles"""
        return latgas_config.config.sum()

    def trialstep(self, latgas_config, energy):
        # energy is just a placeholder and isn't used
        # here
        # choose x,y randomly
        x = rand.randrange(latgas_config.lenX)
        y = rand.randrange(latgas_config.lenY)
        # Position of flipped occupancy to be used by newconfig():
        dconfig = [x, y]
        config = latgas_config.config

        # Calculate energy change if occupation flips at x,y
        x1 = x + 1
        if x == latgas_config.lenX - 1:
            x1 = 0
        y1 = y + 1
        if y == latgas_config.lenY - 1:
            y1 = 0
        dE = -(2 * config[x, y] - 1) * (
            self.J
            * (config[x - 1, y] + config[x1, y] + config[x, y - 1] + config[x, y1])
            + self.Eads
            - self.mu
        )
        return dconfig, dE

    def newconfig(self, latgas_config, dconfig):
        """Construct the new configuration after the trial step is accepted"""
        latgas_config.config[dconfig[0], dconfig[1]] = (
            1 - latgas_config.config[dconfig[0], dconfig[1]]
        )
        return latgas_config


class latgas_config:
    """This class defines the lattice gas model configuration"""

    def __init__(self, lenX, lenY):
        self.lenX = lenX
        self.lenY = lenY
        self.config = np.empty([lenX, lenY])

    def prepare_random(self):
        for i in range(self.lenX):
            for j in range(self.lenY):
                if rand.random() >= 0.5:
                    self.config[i, j] = 1
                else:
                    self.config[i, j] = 0

    def __str__(self):
        s = ""
        for i in range(self.lenX):
            for j in range(self.lenY):
                if self.config[i, j] < 0:
                    s += "-"
                else:
                    s += "+"
            s += "\n"
        return s


class observer(observer_base):
    def __init__(self):
        self.energy_obs = []
        self.density_obs = []

    def logfunc(self, calc_state):
        energy = calc_state.energy
        num = calc_state.model.density(calc_state.config)
        # self.energy_obs.append(energy)
        self.density_obs.append(num)
        return energy, num


if __name__ == "__main__":
    kb = 1.38064852e-23
    m = 1.0 / 6.023 * 1e-26
    T = 300
    h = 6.62607004e-34
    size = 5
    nspin = size * size
    eqsteps = nspin * 1000
    mcsteps = nspin * 1000
    sample_frequency = 1  # nspin
    config = latgas_config(size, size)
    config.prepare_random()
    binning_file = open("binning.dat", "a")
    mu0 = -kb * T * np.log((2.0 * np.pi * m * kb * T / h ** 2) ** (3 / 2) * kb * T)
    Eads = -0.4 * 1.602e-19
    for p in np.linspace(1e4, 0.0001, 20):
        mu = mu0 + kb * T * np.log(p)
        model = latticegas(Eads=Eads, J=0, mu=mu)
        kT = kb * T
        calc = CanonicalMonteCarlo(model, kT, config)
        calc.run(eqsteps)
        myobserver = observer()
        obs = calc.run(mcsteps, sample_frequency, myobserver)
        print(p, "\t", "\t".join([str(x / nspin) for x in obs]))
        # binning analysis
        error_estimate = binning(myobserver.density_obs, 10)
        binning_file.write("\n".join([str(x) for x in error_estimate]) + "\n\n")
        sys.stdout.flush()
        model = calc.model
        config = calc.config
        # print(config)
