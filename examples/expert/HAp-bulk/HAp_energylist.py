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
import sys, os
import copy
from mpi4py import MPI

from pymatgen import Structure
from pymatgen.io.vasp import VaspInput
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen

from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.run_base_mpi import runner

from abics.mc import CanonicalMonteCarlo
from abics.mc_mpi import TemperatureRX_MPI

from model_setup import *


if __name__ == "__main__":
    args = sys.argv
    kB = 8.6173e-5
    nreplicas = int(args[1])
    nprocs_per_vasp = int(args[2])
    commworld = MPI.COMM_WORLD
    worldrank = commworld.Get_rank()
    worldprocs = commworld.Get_size()
    if worldprocs > nreplicas:
        if worldrank == 0:
            print(
                "Setting number of replicas smaller than MPI processes; I hope you"
                + " know what you're doing..."
            )
            sys.stdout.flush()
        if worldrank >= nreplicas:
            # belong to comm that does nothing
            comm = commworld.Split(color=1, key=worldrank)
            comm.Free()
            sys.exit()  # Wait for MPI_finalize
        else:
            comm = commworld.Split(color=0, key=worldrank)
    else:
        comm = commworld

    # prepare config
    cellsize = [2, 2, 2]
    base_structure = Structure.from_file(
        os.path.join(os.path.dirname(__file__), "POSCAR")
    )  # .get_primitive_structure(tolerance=0.001)
    config = HAp_config(base_structure=base_structure, cellsize=cellsize)
    config.set_latgas()

    configs = []
    for i in range(nreplicas):
        configs.append(copy.deepcopy(config))

    # prepare vasp spinel model
    path_to_vasp = "/home/i0009/i000900/src/vasp.5.3/vasp.spawnready.gamma"
    solver = VASPSolver(path_to_vasp)
    vasprun = runner(
        base_input_dir="./baseinput",
        Solver=solver,
        nprocs_per_solver=nprocs_per_vasp,
        comm=MPI.COMM_SELF,
    )

    baseinput = VaspInput.from_directory(
        "baseinput"
    )  # (os.path.join(os.path.dirname(__file__), "baseinput"))
    ltol = 0.1
    #    matcher = StructureMatcher(ltol=ltol, primitive_cell=False, ignored_species=["Pt"])
    #    matcher_site = StructureMatcher(ltol=ltol, primitive_cell=False,stol=0.5,
    #                                             allow_subset=True)#,
    # comparator=FrameworkComparator(), ignored_species=["Pt","Zr"])
    drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
    queen = BorgQueen(drone)
    model = dft_HAp(
        calcode="VASP", vasp_run=vasprun, base_vaspinput=baseinput, queen=queen
    )
    # matcher=matcher, matcher_site=matcher_site, queen=queen, selective_dynamics=["Pt"])

    if worldrank == 0:
        print(config.structure)
        # print(model.xparam(config))

    kTs = kB * np.array([500.0 * 1.1 ** i for i in range(nreplicas)])
    # configs = pickle.load(open("config.pickle","rb"))
    # configs = [copy.deepcopy(config) for i in range(nreplicas)]
    RXcalc = TemperatureRX_MPI(comm, CanonicalMonteCarlo, model, configs, kTs)
    # RXcalc.reload()
    obs = RXcalc.run(
        nsteps=10,
        RXtrial_frequency=2,
        sample_frequency=1,
        observfunc=observables,
        subdirs=True,
    )
    if worldrank == 0:
        for i in range(len(kTs)):
            with open("T" + str(i) + ".dat", "w") as f:
                f.write("\n".join([str(obs[i, j]) for j in range(len(obs[i, :]))]))

    for i in range(9):
        RXcalc.reload()
        obs += RXcalc.run(
            nsteps=10,
            RXtrial_frequency=2,
            sample_frequency=1,
            observfunc=observables,
            subdirs=True,
        )
        obs_write = obs / float(i + 2)
        if worldrank == 0:
            for i in range(len(kTs)):
                with open("T" + str(i) + ".dat", "w") as f:
                    f.write(
                        "\n".join(
                            [str(obs_write[i, j]) for j in range(len(obs_write[i, :]))]
                        )
                    )
