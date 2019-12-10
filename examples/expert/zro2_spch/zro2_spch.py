import numpy as np
import random as rand
import sys, os
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

    # prepare spinel_config
    cellsize = 1
    base_structure = Structure.from_file(
        os.path.join(os.path.dirname(__file__), "POSCAR")
    )  # .get_primitive_structure(tolerance=0.001)
    config = spch_config(base_structure=base_structure, N_Ovac=0.1)
    # config.prepare_Ovac(0.1)
    # config.prepare_ordered()

    configs = []
    for i in range(nreplicas):
        config.prepare_Ovac()
        configs.append(copy.deepcopy(config))

    # prepare vasp spinel model
    vasprun = vasp_run_mpispawn(
        "/home/i0009/i000900/src/vasp.5.3/vasp.spawnready.gamma",
        nprocs=nprocs_per_vasp,
        comm=comm,
    )
    baseinput = VaspInput.from_directory(
        "baseinput"
    )  # (os.path.join(os.path.dirname(__file__), "baseinput"))
    ltol = 0.1
    matcher = StructureMatcher(ltol=ltol, primitive_cell=False, ignored_species=["Pt"])
    matcher_site = StructureMatcher(
        ltol=ltol, primitive_cell=False, stol=0.5, allow_subset=True
    )  # ,
    # comparator=FrameworkComparator(), ignored_species=["Pt","Zr"])
    drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
    queen = BorgQueen(drone)
    model = dft_zro2_spch(
        calcode="VASP",
        vasp_run=vasprun,
        base_vaspinput=baseinput,
        matcher=matcher,
        matcher_site=matcher_site,
        queen=queen,
        selective_dynamics=["Pt"],
    )

    if worldrank == 0:
        print(config.structure)
        # print(model.xparam(config))

    kTs = kB * np.array([500.0 * 1.1 ** i for i in range(nreplicas)])
    # configs = pickle.load(open("config.pickle","rb"))
    # configs = [copy.deepcopy(config) for i in range(nreplicas)]
    RXcalc = TemperatureRX_MPI(comm, CanonicalMonteCarlo, model, configs, kTs)
    RXcalc.reload()
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
