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
import os
import copy
from mpi4py import MPI

from pymatgen import Structure
from abics.mc import CanonicalMonteCarlo, grid_1D, obs_encode, obs_decode
from abics.mc_mpi import RX_MPI_init, TemperatureRX_MPI
from abics.applications.latgas_abinitio_interface.model_setup import (
    group,
    defect_sublattice,
    config,
    dft_latgas
)
from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.run_base_mpi import runner

def g_r(structure, specie1, specie2, grid_1D):
    X = grid_1D.x
    dr = grid_1D.dx

    lattice = structure.lattice
    types_of_specie = [element.symbol for element in structure.types_of_specie]
    assert specie1 in types_of_specie
    assert specie2 in types_of_specie

    structure1 = structure.copy()
    structure2 = structure.copy()

    # print(specie1)
    not_specie1 = copy.copy(types_of_specie)
    not_specie1.remove(specie1)
    not_specie2 = types_of_specie
    not_specie2.remove(specie2)
    # print(not_specie1,not_specie2)

    structure1.remove_species(not_specie1)
    structure2.remove_species(not_specie2)

    num_specie1 = structure1.num_sites
    num_specie2 = structure2.num_sites

    dist = lattice.get_all_distances(structure1.frac_coords, structure2.frac_coords)
    dist_bin = np.around(dist / dr)
    # print(num_specie1,num_specie2)
    g = np.zeros(len(X))
    for i in range(len(X)):
        g[i] = np.count_nonzero(dist_bin == i + 1)

    if specie1 == specie2:
        pref = lattice.volume / (
            4.0 * np.pi * X * X * dr * num_specie1 * (num_specie1 - 1)
        )
    else:
        pref = lattice.volume / (4.0 * np.pi * X * X * dr * num_specie1 * num_specie2)

    g *= pref
    return g

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
path_to_vasp="/home/i0009/i000900/src/vasp.5.3/vasp.spawnready.gamma"
solver = VASPSolver(path_to_vasp)
energy_calculator = runner(
    base_input_dir="./baseinput",
    Solver=solver,
    nprocs_per_solver=nprocs_per_replica,
    comm=MPI.COMM_SELF,
    perturb=0.1
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
