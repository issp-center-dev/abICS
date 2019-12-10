import numpy as np
import random as rand
import sys
import os
import copy
import pickle
from mpi4py import MPI

from pymatgen import Lattice, Structure, Element, PeriodicSite
from pymatgen.io.vasp import Poscar, VaspInput
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
from py_mc.mc import CanonicalMonteCarlo, grid_1D, obs_encode, obs_decode
from py_mc.mc_mpi import RX_MPI_init, TemperatureRX_MPI
from py_mc.applications.latgas_abinitio_interface.model_setup import (
    group,
    defect_sublattice,
    config,
    dft_latgas,
    g_r,
)
from py_mc.applications.latgas_abinitio_interface.run_vasp_mpi import (
    test_runner,
    vasp_runner,
)


def observables(MCcalc, outputfi):
    energy = MCcalc.energy
    nup = MCcalc.config.count("OH", 0)[0]
    ndown = MCcalc.config.count("OH", 1)[0]
    tot_pol = np.abs(nup - ndown) / 6
    grid = grid_1D(0.02, 0.02, 5.0)
    g_r_dat = g_r(MCcalc.config.structure, "O", "H", grid)
    # energy2 = energy**2.0
    # xparam = MCcalc.model.xparam(MCcalc.config)
    outputfi.write(
        "\t".join(
            [str(observable) for observable in [MCcalc.kT, energy, nup, ndown, tot_pol]]
        )
        + "\n"
    )
    outputfi.flush()
    scalar_obs = [MCcalc.kT, energy, nup, ndown, tot_pol]
    return obs_encode(scalar_obs, g_r_dat)


kB = 8.6173e-5
comm, nreplicas, nprocs_per_replica = RX_MPI_init()

# Set up the defect sublattice and construct the HAp configurations

# defect sublattice setup
################## should only need to edit here #################
site_centers = [[0.0, 0.0, 0.25], [0.0, 0.0, 0.75]]
vac = group("V", [], [])
OHgroup = group(
    "OH",
    ["O", "H"],
    np.array(
        [
            [[0.0, 0.0, -0.8438 / 2 + 0.878], [0.0, 0.0, 0.8438 / 2 + 0.878]],
            [[0.0, 0.0, 0.8438 / 2 - 0.878], [0.0, 0.0, -0.8438 / 2 - 0.878]],
        ]
    ),
)
Ogroup = group("O", ["O"])
# groups = [vac, OHgroup, Ogroup]
groups = [OHgroup]
kTstart = 1000.0
kTstep = 1.2
eqsteps = 100
nsteps = 10000
RXtrial_frequency = 2
sample_frequency = 1
dr = 0.01
maxr = 5
#################################################################

#################### config setup ###############################
OH_lattice = defect_sublattice(site_centers, groups)
# num_defects = {"V":1, "OH":4, "O":1}
num_defects = {"OH": 6}
base_structure = Structure.from_file(os.path.join(os.path.dirname(__file__), "POSCAR"))
HAp_config = config(base_structure, OH_lattice, num_defects, [1, 1, 3])
HAp_config.shuffle()
configs = []
for i in range(nreplicas):
    configs.append(copy.deepcopy(HAp_config))
################################################################

# print(HAp_config.structure)
# poscar = Poscar(HAp_config.structure)
# poscar.write_file("POSCAR.vasp")

################### model setup ###############################
# baseinput = VaspInput.from_directory("baseinput")
energy_calculator = vasp_runner(
    base_input_dir="./baseinput",
    path_to_vasp="/home/i0009/i000900/src/vasp.5.3/vasp.spawnready.gamma",
    nprocs_per_vasp=nprocs_per_replica,
    comm=MPI.COMM_SELF,
    perturb=0.1,
)
# energy_calculator = test_runner()
model = dft_latgas(energy_calculator, save_history=True)
##############################################################

################### RXMC calculation #########################
kTs = kB * np.array([kTstart * kTstep ** i for i in range(nreplicas)])
grid = grid_1D(dr, dr, maxr)
RXcalc = TemperatureRX_MPI(comm, CanonicalMonteCarlo, model, configs, kTs)

obs = RXcalc.run(
    eqsteps, RXtrial_frequency, sample_frequency, observfunc=observables, subdirs=True
)
# RXcalc.reload()
obs = RXcalc.run(
    nsteps, RXtrial_frequency, sample_frequency, observfunc=observables, subdirs=True
)


if comm.Get_rank() == 0:
    for i in range(len(kTs)):
        scalar_obs, g_OH = obs_decode(obs[i])
        grid = grid_1D(0.02, 0.02, 5.0)
        g_r_r = zip(grid.x, g_OH)
        with open("grT" + str(i) + ".dat", "w") as f:
            for dat in g_r_r:
                f.write("\t".join([str(xy) for xy in dat]) + "\n")
        print("\t".join([str(x) for x in scalar_obs]))
