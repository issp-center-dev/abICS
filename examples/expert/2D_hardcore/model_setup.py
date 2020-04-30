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
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from abics.mc import model

def gauss(x, x0, sigma):
    return (
        1.0
        / (np.sqrt(2.0 * np.pi) * sigma)
        * np.exp(-np.power((x - x0) / sigma, 2.0) / 2.0)
    )


class HC_2D(model):
    """This class defines the 2D hardcore potential model"""

    model_name = "HC_2D"

    def __init__(self, maxstep):
        self.maxstep = maxstep

    def energy(self, HC_2D_config):
        """ Return a very large number if there is any overlap, else 0"""
        coors = HC_2D_config.coors
        Lcell = HC_2D_config.Lcell
        # make use of newaxis and broadcasting to make NxN array of
        # distance vectors sij[N,N,ndim]:
        sij = coors[:, None] - coors[None]
        sij = sij - np.around(sij / Lcell) * Lcell  # minimum image convention
        dist = np.linalg.norm(sij, axis=2)
        dist = squareform(dist) - HC_2D_config.rdisc * 2  # condensed format
        if np.all(np.greater(dist, 0)):
            return 0
        else:
            return float("Infinity")

    def trialstep(self, HC_2D_config, energy_now):

        Lcell = HC_2D_config.Lcell
        rdisc = HC_2D_config.rdisc
        # Randomly choose a disc
        move_id = np.random.randint(HC_2D_config.Ndisc)
        coors = HC_2D_config.coors
        move_coor = coors[move_id].copy()
        move_coor += self.maxstep * np.array(rand.uniform(-1, 1))
        trial_coor_others = np.delete(coors, move_id, 0)

        dconfig = [move_id, move_coor]

        sij = trial_coor_others - move_coor
        sij = sij - np.around(sij / Lcell) * Lcell
        dist = np.linalg.norm(sij, axis=1) - rdisc * 2

        if np.all(np.greater(dist, 0)):
            dE = 0
        else:
            dE = float("Infinity")

        return dconfig, dE

    def newconfig(self, HC_2D_config, dconfig):
        """Construct the new configuration after the trial step is accepted"""
        HC_2D_config.coors[dconfig[0]] = dconfig[1]
        return HC_2D_config


class HC_2D_config:
    """This class defines the 2D hardcore discs config"""

    def __init__(self, Lcell, Ndisc, rdisc):
        self.Lcell = Lcell
        self.Ndisc = Ndisc
        self.rdisc = rdisc
        self.coors = np.zeros((Ndisc, 2))
        self.g_pref = Lcell ** 2.0 / (2.0 * np.pi * Ndisc * (Ndisc - 1))

    def prepare_ordered(self):
        # Prepare triangular lattice with a little space
        # First filled along Y direction (will be space in X dir)

        disc_spacing_x = 1.05 * self.rdisc * 2.0
        disc_spacing_y = disc_spacing_x * 3.0 ** 0.5 * 0.5
        maxX = int(self.Lcell / disc_spacing_x)
        maxY = int(self.Lcell / disc_spacing_y)
        maxX2, remainder = divmod(self.Ndisc, maxY)
        assert maxX2 < maxX

        id = 0
        for i in range(maxX2):
            for j in range(maxY):
                self.coors[id, 0] = i + 0.5 * j % 1
                self.coors[id, 1] = j
                id += 1
        for i in range(remainder):
            self.coors[id, 0] = maxX2 + 0.5 * i % 1
            self.coors[id, 1] = i
            id += 1
        assert id == self.Ndisc

        self.coors[:, 0] *= disc_spacing_x
        self.coors[:, 1] *= disc_spacing_y


def g_r(HC_2D_config, grid_1D):
    X = grid_1D.x
    dr = grid_1D.dx

    # First construct distance matrix
    coors = HC_2D_config.coors
    Lcell = HC_2D_config.Lcell
    Ndisc = HC_2D_config.Ndisc
    # make use of newaxis and broadcasting to make NxN array of
    # distance vectors sij[N,N,ndim]:
    sij = coors[:, None] - coors[None]
    sij = sij - np.around(sij / Lcell) * Lcell  # minimum image convention
    dist = squareform(np.linalg.norm(sij, axis=2))
    dist_bin = np.around(dist / dr)
    g = np.zeros(len(X))
    for i in range(len(X)):
        g[i] = np.count_nonzero(dist_bin == i + 1) * 2

    g *= HC_2D_config.g_pref / (X * dr)
    return g


def observables(MCcalc, outputfi):
    return g_r(MCcalc.config, MCcalc.grid)


def plot_fig(MCcalc, nstep, sample_frequency=100):
    plt.axis([0, MCcalc.config.Lcell, 0, MCcalc.config.Lcell])
    plt.ion()
    for i in range(nstep):
        MCcalc.MCstep()
        if i % sample_frequency == 0:
            for i in range(MCcalc.config.Ndisc):
                circle = plt.Circle(
                    MCcalc.config.coors[i], radius=MCcalc.config.rdisc, fc="b"
                )
                plt.gca().add_patch(circle)
            plt.pause(0.1)
