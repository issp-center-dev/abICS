import numpy as np
import random as rand
import sys
import os
import copy
import pickle
from mpi4py import MPI

from pymatgen import Lattice, Structure, Element
from pymatgen.io.vasp import Poscar, VaspInput
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen

# from mc.applications.dft_spinel_mix.dft_spinel_mix import dft_spinel_mix, spinel_config
from applications.dft_spinel_mix.run_vasp_mpi import vasp_run_mpispawn
from mc import (
    model,
    CanonicalMonteCarlo,
    MultiProcessReplicaRun,
    TemperatureReplicaExchange,
)
from mc_mpi import TemperatureRX_MPI

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
    vasprun = vasp_run_mpispawn(
        "/home/i0009/i000900/src/vasp.5.3/vasp.spawnready.gamma",
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
