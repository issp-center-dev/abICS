import numpy as np
import random as rand
import sys, os
import copy
import pickle

# from mpi4py import MPI

from pymatgen import Lattice, Structure, Element
from pymatgen.io.vasp import Poscar, VaspInput
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen

# from mc.applications.dft_spinel_mix.dft_spinel_mix import dft_spinel_mix, spinel_config
# from applications.dft_spinel_mix.run_vasp_mpi import vasp_run_mpispawn
from mc import (
    model,
    CanonicalMonteCarlo,
    MultiProcessReplicaRun,
    TemperatureReplicaExchange,
)

# from mc_mpi import TemperatureRX_MPI

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
