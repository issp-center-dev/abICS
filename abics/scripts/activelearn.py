# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2020- The University of Tokyo
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

from typing import MutableMapping

import sys
import os
import glob
import datetime

import numpy as np

from abics import __version__

from abics.mc_mpi import (
    RX_MPI_init,
    RefParams,
)
from abics.applications.latgas_abinitio_interface import map2perflat
# from abics.applications.latgas_abinitio_interface import DefaultObserver
from abics.applications.latgas_abinitio_interface.model_setup import (
    perturb_structure,
)
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)
from abics.applications.latgas_abinitio_interface.base_solver import SolverBase
from abics.applications.latgas_abinitio_interface.params import ALParams, DFTParams

from abics.util import exists_on_all_nodes
from pymatgen.core import Structure

import logging
import abics.loggers as loggers
from pathlib import Path

#loggers.set_log_handles(app_name="activelearn", level=logging.INFO, mpi_log="master")
loggers.set_log_handles(app_name="activelearn", level=logging.INFO, mpi_log="collect", log_path=Path("run.log"))
logger = logging.getLogger("main")


def main_impl(params_root: MutableMapping):
    rxparams = RefParams.from_dict(params_root["mlref"])
    # nprocs_per_replica = rxparams.nprocs_per_replica
    # nreplicas = rxparams.nreplicas
    # nsteps = rxparams.nsteps
    # sample_frequency = rxparams.sample_frequency
    ndata = rxparams.ndata
    comm = RX_MPI_init(rxparams.nreplicas, rxparams.seed)
    alparams = ALParams.from_dict(params_root["mlref"]["solver"])
    mcparams = DFTParams.from_dict(params_root["sampling"]["solver"])
    configparams = DFTConfigParams.from_dict(params_root["config"])

    solver: SolverBase = SolverBase.create(alparams.solver, alparams)

    perturb = alparams.perturb
    solver_output = solver.output
    solver_input = solver.input

    myreplica = comm.Get_rank()

    # Find newest MC run
    nextMC_index = 0
    while exists_on_all_nodes(comm, "MC{}".format(nextMC_index)):
        nextMC_index += 1

    rootdir = os.getcwd()
    if nextMC_index == 0:  # "MC0" does not exist
        # Random sampling!
        ndigits = len(str(ndata))
        fmtstr = "input{:0>" + str(ndigits) + "d}"

        if exists_on_all_nodes(comm, "AL0"):
            if exists_on_all_nodes(comm, "ALloop.progress"):
                logger.error("It seems you've already run the first active learning step. You should train now.")
                sys.exit(1)
            os.chdir("AL0")

            # attempting post processing of ALed output
            if not exists_on_all_nodes(comm, "baseinput.progress"):  # 1st step
                runstep = 0

            else:
                with open("baseinput.progress", "r") as fi:
                    runstep = int(fi.readlines()[-1]) + 1

            if runstep + 1 == len(alparams.base_input_dir):
                # we're done with this AL step
                finalrun = True
            else:
                finalrun = False
                solver_input.from_directory(alparams.base_input_dir[runstep + 1])
            energies = []

            # Preparing ignored species structure to add in postprocess
            config = defect_config(configparams)
            if alparams.ignore_species:
                ignore_structure = config.dummy_structure()
                remove_sp = filter(
                    lambda sp: sp not in alparams.ignore_species,
                    ignore_structure.symbol_set,
                )
                ignore_structure.remove_species(remove_sp)
            rundir_list = []

            logger.info(f"-Parsing {alparams.solver} results in AL0/*/input*/baseinput{runstep}...")
            if finalrun:
                logger.info(f"--This is the final {alparams.solver} calculation step for AL0")
            else:
                logger.info("--Input files for the next calculation step will be")
                logger.info(f"---prepared in AL{nextMC_index}/*/input*/baseinput{runstep+1}")

            for i in range(ndata):
                energy, st = solver_output.get_results(
                    os.path.join(
                        str(myreplica), fmtstr.format(i), f"baseinput{runstep}"
                    )
                )
                if finalrun:
                    # energy_calculator may return the structure w/o ignored structure
                    if alparams.ignore_species:
                        for site in ignore_structure:
                            st.append(site.species_string, site.frac_coords)
                    st.sort(key=lambda site: site.species_string)
                    st.to(
                        fmt="POSCAR",
                        filename=os.path.join(str(myreplica), f"structure.{i}.vasp"),
                    )
                    energies.append([energy])
                else:
                    solver_input.update_info_by_structure(st)
                    inputfile = os.path.abspath(
                        os.path.join(
                            str(myreplica),
                            fmtstr.format(i),
                            f"baseinput{runstep+1}",
                        )
                    )
                    solver_input.write_input(inputfile)
                    rundir_list.append(inputfile)
            rundir_list = comm.gather(rundir_list, root=0)
            if myreplica == 0 and not finalrun:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                with open(os.path.join(rootdir, "rundirs.txt"), "w") as fi:
                    fi.write("\n".join(rundir_list))
                    fi.write("\n")
                    fi.flush()
                    os.fsync(fi.fileno())
            if finalrun:
                with open(os.path.join(str(myreplica), "energy_corr.dat"), "w") as fi:
                    for i in range(len(energies)):
                        energy = energies[i][0]
                        fi.write(f"NaN {energy} {i}\n")
                os.chdir(rootdir)
                comm.Barrier()
                if myreplica == 0:
                    with open("ALloop.progress", "a") as fi:
                        logger.info("-Writing ALloop.progress")
                        fi.write("AL0\n")
                        fi.flush()
                        os.fsync(fi.fileno())
                    now = datetime.datetime.now()
                    logger.info("--Done. Please run abics_train next.")
                    logger.info(f"Exiting normally on {now}.\n")
                sys.exit(0)
            else:
                comm.Barrier()
                if myreplica == 0:
                    with open("baseinput.progress", "a") as fi:
                        fi.write(str(runstep) + "\n")
                        fi.flush()
                        os.fsync(fi.fileno())

                    logger.info(f"-Finished preparing {alparams.solver} input in AL0/*/input*/baseinput{runstep+1}.")
                    logger.info("--See rundirs.txt for a full list of the directories.")
                    logger.info("-Please perform the calculations in those directories before running abics_mlref again.")
                    now = datetime.datetime.now()
                    logger.info(f"Exiting normally on {now}.\n")
                sys.exit(0)
        else:
            # "AL0" nor "MC0" does not exist:
            # this should be the first time running this script
            if myreplica == 0:
                logger.info("-No preceding MC or AL runs exist. Preparing {} input by random sampling".format(alparams.solver))
                logger.info("--Input files will be prepared in AL0/*/input*/baseinput0")

            configparams = DFTConfigParams.from_dict(params_root["config"])
            config = defect_config(configparams)
            comm.Barrier()
            if myreplica == 0:
                os.mkdir("AL0")

            # Wait until the new directory is visible from all ranks
            if not exists_on_all_nodes(comm, "AL0"):
                logger.error("Failed to create AL0 directory.")
                sys.exit(1)

            os.chdir("AL0")
            rundir_list = []
            solver_input = solver.input
            solver_input.from_directory(alparams.base_input_dir[0])
            for i in range(ndata):
                config.shuffle() # randomize config
                config.structure.sort(key=lambda site: site.species_string)
                structure0 = config.structure.copy()  
                perturb_structure(config.structure, perturb)
                solver_input.update_info_by_structure(config.structure)
                inputfile = os.path.abspath(
                    os.path.join(str(myreplica), fmtstr.format(i), "baseinput0")
                )
                solver_input.write_input(inputfile)
                rundir_list.append(inputfile)

                # save structure before running through solver (contains ignored species)
                structure0.to(
                    fmt='POSCAR',
                    filename=os.path.abspath(
                        os.path.join(str(myreplica), fmtstr.format(i), "initial.vasp")
                        ),
                )
            rundir_list = comm.gather(rundir_list, root=0)
            if myreplica == 0:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                with open(os.path.join(rootdir, "rundirs.txt"), "w") as fi:
                    fi.write("\n".join(rundir_list))
                    fi.write("\n")
                    fi.flush()
                    os.fsync(fi.fileno())

                logger.info(f"-Finished preparing {alparams.solver} input in AL0/*/input*/baseinput0.")
                logger.info("--See rundirs.txt for a full list of the directories.")
                logger.info("--Please perform the calculations in those directories before running abics_mlref again.")
                now = datetime.datetime.now()
                logger.info(f"Exiting normally on {now}.\n")

                sys.exit(0)

    else:  # "MC*" exists
        # Active learning!
        MCdirname = f"MC{nextMC_index-1}"
        MCdir = os.path.join(os.getcwd(), MCdirname)
        nsamples = 0
        with open(os.path.join(MCdir, str(myreplica), "obs.dat")) as f:
            for _ in f:
                nsamples += 1
        ndigits = len(str(nsamples))
        fmtstr = "input{:0>" + str(ndigits) + "d}"
        with open("ALloop.progress", "r") as fi:
            last_li = fi.readlines()[-1].strip()
        if last_li != MCdirname:
            logger.error("The last action:")
            logger.error(f"  expected: {MCdirname}")
            logger.error(f"  actual:   {last_li}")
            logger.error("You shouldn't run activelearn now.")
            logger.error("Either train (abics_train) or MC (abics_sampling) first.")
            sys.exit(1)
        obs = np.load(os.path.join(MCdir, str(myreplica), "obs_save.npy"))
        energy_ref = obs[:, 0]
        ALstep = nextMC_index
        ALdir = os.path.join(os.getcwd(), f"AL{ALstep}", str(myreplica))
        config = defect_config(configparams)

        if exists_on_all_nodes(comm, ALdir):  # We're going to read solver results
            os.chdir(ALdir)

            if not exists_on_all_nodes(
                comm, os.path.join(os.path.pardir, "baseinput.progress")
            ):
                # 1st runner step
                runstep = 0
            else:
                # 2nd or later runner step
                with open(os.path.join(os.pardir, "baseinput.progress"), "r") as fi:
                    runstep = int(fi.readlines()[-1]) + 1

            if runstep + 1 == len(alparams.base_input_dir):
                # This is the last run for this AL step
                finalrun = True

            else:
                finalrun = False
                solver_input.from_directory(alparams.base_input_dir[runstep + 1])

            logger.info(f"-Parsing {alparams.solver} results in AL{nextMC_index}/*/input*/baseinput{runstep}...")
            if finalrun:
                logger.info(f"--This is the final {alparams.solver} calculation step for AL{nextMC_index}")
            else:
                logger.info("--Input files for the next calculation step will be")
                logger.info(f"---prepared in AL{nextMC_index}/*/input*/baseinput{runstep+1}")

            energy_corrlist = []
            relax_max = []
            rundir_list = []
            sample_list = []
            for path in glob.iglob(os.path.join(ALdir, "input*")):
                sample_list.append(int(os.path.basename(path).split("t")[1]))
            sample_list.sort()
            for i in sample_list:
                energy, st_rel = solver_output.get_results(
                    os.path.join(fmtstr.format(i), f"baseinput{runstep}")
                )
                if finalrun:
                    # Get original structure for calculating relaxation magnitude
                    st_in = Structure.from_file(
                        os.path.join(fmtstr.format(i), "initial.vasp")
                    )
                    
                    # energy_calculator may return the structure w/o ignored structure
                    # so we need to add them
                    if alparams.ignore_species:
                        # initial.vasp contains full structure with ignored atoms
                        ignore_structure = st_in.copy()
                        remove_sp = filter(
                            lambda sp: sp not in alparams.ignore_species, ignore_structure.symbol_set
                            )
                        ignore_structure.remove_species(remove_sp)
                        for site in ignore_structure:
                            st_rel.append(site.species_string, site.frac_coords)

                    st_rel.to(fmt="POSCAR", filename=f"structure.{i}.vasp")
                    energy_corrlist.append([energy_ref[i], energy, i])

                    # Examine relaxation
                    dmax = 0
                    dmax_id = 0
                    for j, site in enumerate(st_in):
                        d = site.distance(st_rel[j])
                        if d > dmax:
                            dmax = d
                            dmax_id = j
                    relax_max.append([dmax_id, dmax])

                else:
                    solver_input.update_info_by_structure(st_rel)
                    filename = os.path.abspath(
                        os.path.join(fmtstr.format(i), f"baseinput{runstep+1}")
                    )
                    solver_input.write_input(filename)
                    rundir_list.append(filename)
            rundir_list = comm.gather(rundir_list, root=0)
            if myreplica == 0 and not finalrun:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                with open(os.path.join(rootdir, "rundirs.txt"), "w") as fi:
                    fi.write("\n".join(rundir_list))
                    fi.write("\n")
                    fi.flush()
                    os.fsync(fi.fileno())

            if finalrun:
                with open("energy_corr.dat", "w") as fi:
                    for i in range(len(energy_corrlist)):
                        ene_nn, ene_ref, step = energy_corrlist[i]
                        fi.write(f"{ene_nn} {ene_ref} {step}\n")
                with open("relax_max.dat", "w") as fi:
                    for row in relax_max:
                        fi.write("{} \t {}\n".format(row[0], row[1]))
                os.chdir(rootdir)
                comm.Barrier()
                if myreplica == 0:
                    with open("ALloop.progress", "a") as fi:
                        logger.info("-Writing ALloop.progress")
                        fi.write("AL{}\n".format(ALstep))
                        fi.flush()
                        os.fsync(fi.fileno())
                    now = datetime.datetime.now()
                    logger.info("--Done. Please run abics_train next.")
                    logger.info(f"Exiting normally on {now}.\n")
                sys.exit(0)
            else:
                comm.Barrier()
                if myreplica == 0:
                    with open(
                        os.path.join(os.path.pardir, "baseinput.progress"), "a"
                    ) as fi:
                        fi.write(str(runstep) + "\n")
                        fi.flush()
                        os.fsync(fi.fileno())
                    logger.info(f"-Finished preparing {alparams.solver} input in AL{nextMC_index}/*/input*/baseinput{runstep+1}.")
                    logger.info("--See rundirs.txt for a full list of the directories.")
                    logger.info("--Please perform the calculations in those directories before running abics_mlref again.")
                    now = datetime.datetime.now()
                    logger.info(f"Exiting normally on {now}.\n")
                sys.exit(0)

        else:  # Create first input for this AL step
            logger.info(f"-This is the first run for AL{nextMC_index}")
            logger.info(f"--Configurations from MC{nextMC_index-1} will be converted to {alparams.solver} input")
            os.makedirs(ALdir, exist_ok=False)
            solver_input.from_directory(alparams.base_input_dir[0])
            os.chdir(ALdir)
            sample_list = rxparams.sampling(nsamples)
            rundir_list = []
            for i in sample_list:
                # We read structures from MC, but we also need to add ignored species.
                # We utilize structure_norel.*.vasp, which is the structure before solver run
                # in MC with ignored species.
                st = Structure.from_file(
                    os.path.join(MCdir, str(myreplica), f"structure.{i}.vasp")
                )
                ignore_structure = Structure.from_file(
                    os.path.join(MCdir, str(myreplica), f"structure_norel.{i}.vasp")
                )

                if mcparams.ignore_species:
                    remove_sp = filter(
                        lambda sp: sp not in mcparams.ignore_species, ignore_structure.symbol_set
                        )
                    ignore_structure.remove_species(remove_sp)
                    for site in ignore_structure:
                        st.append(site.species_string, site.frac_coords)
                st.sort(key=lambda site: site.species_string)
                st0 = st.copy()
                
                perturb_structure(st, perturb)
                solver_input.update_info_by_structure(st)
                solver_input.write_input(
                    os.path.abspath(os.path.join(fmtstr.format(i), "baseinput0"))
                )
                st0.to(fmt='POSCAR',filename=os.path.abspath(os.path.join(fmtstr.format(i), "initial.vasp")))
                rundir_list.append(
                    os.path.abspath(os.path.join(fmtstr.format(i), "baseinput0"))
                )
            rundir_list = comm.gather(rundir_list, root=0)
            if myreplica == 0:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                with open(os.path.join(rootdir, "rundirs.txt"), "w") as fi:
                    fi.write("\n".join(rundir_list))
                    fi.write("\n")
                    fi.flush()
                    os.fsync(fi.fileno())
                logger.info(f"-Finished preparing {alparams.solver} input in AL{nextMC_index}/*/input*/baseinput0.")
                logger.info("--See rundirs.txt for a full list of the directories.")
                logger.info("--Please perform the calculations in those directories before running abics_mlref again.")
                now = datetime.datetime.now()
                logger.info(f"Exiting normally on {now}.\n")
            sys.exit(0)


def main():
    import toml

    now = datetime.datetime.now()
    logger.info(f"Running abics_mlref (abICS v{__version__}) on {now}")

    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    logger.info(f"-Reading input from: {tomlfile}")
    params_root = toml.load(tomlfile)
    main_impl(params_root)


if __name__ == "__main__":
    main()
