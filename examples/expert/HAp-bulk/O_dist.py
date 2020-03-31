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
import copy
import pickle
from model_setup import *


def gauss(x, x0, sigma):
    return (
        1.0
        / (np.sqrt(2.0 * np.pi) * sigma)
        * np.exp(-np.power((x - x0) / sigma, 2.0) / 2.0)
    )


if __name__ == "__main__":
    X, Y = np.mgrid[0:23:0.05, 0:1]
    sigma = 1.0
    data = np.vstack((X.flatten(), Y.flatten())).T
    data_base = copy.deepcopy(data)
    #    T_to_rank = pickle.load
    config = pickle.load(open(str(6) + "/calc.pickle", "rb"))
    config = copy.deepcopy(config)
    structure = config.structure
    structure = config.base_structure
    structure.remove_species(["Pt", "Zr"])
    # structure_base.remove_species(["Pt","Zr"])
    for site in structure:
        x0 = site.coords[2]
        for grid in data:
            grid[1] += gauss(grid[0], x0, sigma)
    #    for site in structure_base:
    #        x0 = site.coords[2]
    #        for grid in data:
    #            grid[1] -= gauss(grid[0], x0, sigma)

    for grid in data:
        print(grid[0], grid[1])
