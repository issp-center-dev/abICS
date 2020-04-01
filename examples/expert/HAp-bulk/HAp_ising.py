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
import sys
import os
import pickle
from mpi4py import MPI

from pymatgen import Structure
from pymatgen.io.vasp import VaspInput
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.run_base_mpi import runner

from model_setup import *
from make_ising_list import *

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
    cellsize = [1, 1, 3]
    base_structure = Structure.from_file(
        os.path.join(os.path.dirname(__file__), "POSCAR")
    )  # .get_primitive_structure(tolerance=0.001)
    config = HAp_config(base_structure=base_structure, cellsize=cellsize)

    # prepare lattice-gas type representation
    reps = make_latgas_reps(spec_id=[0, 1, 2], nspec=[1, 4, 1], spec_orient=[1, 2, 1])
    energy_reps = np.zeros(len(reps))

    # prepare vasp model
    path_to_vasp = "/home/i0009/i000900/src/vasp.5.3/vasp.spawnready.gamma"
    solver = VASPSolver(path_to_vasp)
    vasprun = runner(
        Solver=solver,
        nprocs=nprocs_per_vasp,
        comm=MPI.COMM_SELF,
    )
    baseinput = VaspInput.from_directory(
        "baseinput"
    )  # (os.path.join(os.path.dirname(__file__), "baseinput"))
    ltol = 0.1
    matcher = StructureMatcher(ltol=ltol, primitive_cell=False)
    matcher_base = StructureMatcher(
        ltol=ltol, primitive_cell=False, stol=0.5, allow_subset=True
    )  # ,
    # comparator=FrameworkComparator(), ignored_species=["Pt","Zr"])
    drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
    queen = BorgQueen(drone)
    model = dft_HAp(
        calcode="VASP",
        vasp_run=vasprun,
        base_vaspinput=baseinput,
        matcher_base=matcher_base,
        queen=queen,
        matcher=matcher,
    )
    # matcher=matcher, matcher_site=matcher_site, queen=queen, selective_dynamics=["Pt"])

    myrank = comm.Get_rank()
    os.mkdir(str(myrank))
    os.chdir(str(myrank))
    for i in range(len(reps)):
        if i % nreplicas == myrank:
            config.set_latgas(reps[i])
            print(reps[i])
            energy_reps[i] = model.energy(config, save_history=True)

    energy_buffer = np.empty(len(reps))
    comm.Allreduce(energy_reps, energy_buffer, op=MPI.SUM)
    os.chdir("../")
    if myrank == 0:
        with open("energy_reps.pickle", "wb") as f:
            pickle.dump(energy_buffer, f)

        with open("latgas_reps.pickle", "wb") as f:
            pickle.dump(reps, f)

    print(energy_buffer)
