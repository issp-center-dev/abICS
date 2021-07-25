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

import sys
import os,copy

from mpi4py import MPI
import numpy as np

from abics.mc import RandomSampling

from abics.mc_mpi import (
    RX_MPI_init,
    RefParams,
    EmbarrassinglyParallelSampling,
)
from abics.applications.latgas_abinitio_interface import map2perflat
from abics.applications.latgas_abinitio_interface import default_observer
from abics.applications.latgas_abinitio_interface.model_setup import (
    dft_latgas,
    ObserverParams,
    perturb_structure,
)
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)
from abics.applications.latgas_abinitio_interface.run_base_mpi import (
    runner,
    runner_multistep,
)
from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.qe import QESolver
from abics.applications.latgas_abinitio_interface.aenet import aenetSolver
from abics.applications.latgas_abinitio_interface.openmx import OpenMXSolver
from abics.applications.latgas_abinitio_interface.mocksolver import MockSolver
from abics.applications.latgas_abinitio_interface.params import ALParams

from abics.util import exists_on_all_nodes
from pymatgen import Structure


def main_impl(tomlfile):
    rxparams = RefParams.from_toml(tomlfile)
    nprocs_per_replica = rxparams.nprocs_per_replica
    nreplicas = rxparams.nreplicas
    nsteps = rxparams.nsteps
    sample_frequency = rxparams.sample_frequency
    comm = RX_MPI_init(rxparams)
    alparams = ALParams.from_toml(tomlfile)
    configparams = DFTConfigParams.from_toml(tomlfile)

    if alparams.solver == "vasp":
        solver = VASPSolver(alparams.path, alparams.ignore_species)
    elif alparams.solver == "qe":
        parallel_level = alparams.properties.get("parallel_level", {})
        solver = QESolver(alparams.path, parallel_level=parallel_level)
    elif alparams.solver == "aenet":
        solver = aenetSolver(
            alparams.path, alparams.ignore_species, alparams.solver_run_scheme
        )
    elif alparams.solver == "openmx":
        solver = OpenMXSolver(alparams.path)
    elif alparams.solver == "mock":
        solver = MockSolver()
    else:
        print("unknown solver: {}".format(alparams.solver))
        sys.exit(1)

    perturb = alparams.perturb
    solver_output = solver.output
    solver_input = solver.input
    
    myreplica = comm.Get_rank()
    
    # Find newest MC run
    i = 0
    while exists_on_all_nodes(comm, "MC{}".format(i)):
        i += 1

    rootdir = os.getcwd()
    if i == 0:  # Random sampling!
        ndigits = len(str(nsteps//sample_frequency))
        fmtstr = "input{:0>"+str(ndigits)+"d}"

        if exists_on_all_nodes(comm, "AL0"):
            if exists_on_all_nodes(comm, "ALloop.progress"):
                print(
                    "It seems you've already run the first active learning step. You should train now."
                )
                sys.exit(1)
            os.chdir("AL0")

            # attempting post processing of ALed output
            if not exists_on_all_nodes(comm, "baseinput.progress"): # 1st step
                runstep = 0
            else:
                with open("baseinput.progress", "r") as fi:
                    runstep = int(fi.readlines()[-1]) + 1

            if runstep + 1 == len(alparams.base_input_dir):
                # we're done with this AL step
                finalrun = True
            else:
                finalrun = False
                solver_input.from_directory(alparams.base_input_dir[runstep+1])
            energies = []

            # Preparing ignored species structure to add in postprocess
            config = defect_config(configparams)
            if alparams.ignore_species:
                ignore_structure = config.dummy_structure()
                remove_sp = filter(
                    lambda sp: sp not in alparams.ignore_species, ignore_structure.symbol_set
                )
                ignore_structure.remove_species(remove_sp)
            rundir_list = []
            for i in range(nsteps // sample_frequency):
                energy, st = solver_output.get_results(os.path.join(str(myreplica),fmtstr.format(i),"baseinput{}".format(runstep)))
                if finalrun:
                    # energy_calculator may return the structure w/o ignored structure
                    if alparams.ignore_species:
                        for site in ignore_structure:
                            st.append(site.species_string, site.frac_coords)

                    st.to("POSCAR", os.path.join(str(myreplica),"structure.{}.vasp".format(i)))
                    energies.append([energy])
                else:
                    solver_input.update_info_by_structure(st)
                    solver_input.write_input(
                        os.path.join(str(myreplica),fmtstr.format(i),"baseinput{}".format(runstep+1))
                    )
                    rundir_list.append(
                        os.path.abspath(
                            os.path.join(str(myreplica),fmtstr.format(i),"baseinput{}".format(runstep+1))
                        )
                    )
            rundir_list = comm.gather(rundir_list, root = 0)
            if myreplica == 0 and not finalrun:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                with open(os.path.join(rootdir, "rundirs.txt"), "w") as fi:
                    fi.write("\n".join(rundir_list))
                    fi.flush()
                    os.fsync(fi.fileno())
            if finalrun:
                np.save(os.path.join(str(myreplica),"obs_save.npy"), energies)
                os.chdir(rootdir)
                comm.Barrier()
                if myreplica == 0:
                    with open("ALloop.progress", "a") as fi:
                        fi.write("AL0\n")
                        fi.flush()
                        os.fsync(fi.fileno())
                sys.exit(0)
            else:
                comm.Barrier()
                if myreplica == 0:
                    with open("baseinput.progress", "a") as fi:
                        fi.write(str(runstep) + "\n")
                        fi.flush()
                        os.fsync(fi.fileno())
                sys.exit(0)

            print(
                "It seems you've already run the first active learning step." + \
                "Check if it has completed normally. If it has, then you should train now."
            )
            sys.exit(1)

        # If this line is reached, this should be the first time running this script
        else:
            configparams = DFTConfigParams.from_toml(tomlfile)
            config = defect_config(configparams)
            comm.Barrier()
            if myreplica == 0:
                os.mkdir("AL0")

            # Wait until the new directory is visible from all ranks
            if not exists_on_all_nodes(comm, "AL0"):
                print("Directory creation failed")
                sys.exit(1)

            os.chdir("AL0")
            rundir_list = []
            solver_input = solver.input
            solver_input.from_directory(alparams.base_input_dir[0])
            for i in range(nsteps // sample_frequency):
                config.shuffle()
                perturb_structure(config.structure,perturb)
                config.structure.sort(key=lambda site: site.species_string)
                solver_input.update_info_by_structure(config.structure)
                solver_input.write_input(os.path.join(str(myreplica),fmtstr.format(i),"baseinput0"))
                rundir_list.append(os.path.abspath(os.path.join(str(myreplica),fmtstr.format(i),"baseinput0")))
            rundir_list = comm.gather(rundir_list, root = 0)
            if myreplica == 0:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                with open(os.path.join(rootdir, "rundirs.txt"), "w") as fi:
                    fi.write("\n".join(rundir_list))
                    fi.flush()
                    os.fsync(fi.fileno())

    else:  # Active learning!
        ndigits = len(str(nsteps))
        fmtstr = "input{:0>"+str(ndigits)+"d}"
        MCdir = os.path.join(os.getcwd(), "MC{}".format(i - 1))
        with open("ALloop.progress", "r") as fi:
            last_li = fi.readlines()[-1]
        if "MC{}".format(i - 1) not in last_li:
            print("You shouldn't run activelearn now. Either train or MC first.")
            sys.exit(1)
        MCdir = os.path.join(os.getcwd(), "MC{}".format(i - 1))
        obs = np.load(os.path.join(MCdir, str(myreplica), "obs_save.npy"))
        energy_ref = obs[:, 0]
        ALdir = os.path.join(os.getcwd(), "AL{}".format(i), str(myreplica))
        ALstep = i
        config = defect_config(configparams)
        perf_st = config.dummy_structure()
        if alparams.ignore_species:
            ignore_structure = perf_st.copy()
            remove_sp = filter(
                lambda sp: sp not in alparams.ignore_species, perf_st.symbol_set
            )
            ignore_structure.remove_species(remove_sp)

        
        if exists_on_all_nodes(comm, ALdir): # We're going to read solver results
            os.chdir(ALdir)

            if not exists_on_all_nodes(comm, os.path.join(os.path.pardir, "baseinput.progress")):
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
                solver_input.from_directory(alparams.base_input_dir[runstep+1])
                
            energy_corrlist = []
            relax_max = []
            rundir_list = []
            for i in range(0, nsteps,  sample_frequency):
                energy, st_rel = solver_output.get_results(
                    os.path.join(fmtstr.format(i), "baseinput{}".format(runstep))
                )
                if finalrun:
                    # Get original structure for calculating relaxation magnitude
                    # In st_in, used for aenet,
                    # (i)  "Ignore species" are omitted and should be added.
                    # (ii) "Vacancy element" is included and should be removed.
                    st_in = Structure.from_file(
                        os.path.join(MCdir, str(myreplica), "structure.{}.vasp".format(i))
                    )
                
                    # (i)
                    st = map2perflat(perf_st, st_in)
                    st.remove_species(["X"])
                    # (ii)
                    if alparams.vac_space_holder:
                        st.remove_species(alparams.vac_space_holder)
                    stbak = st.copy()

                    # energy_calculator may return the structure w/o ignored structure
                    if alparams.ignore_species:
                        for site in ignore_structure:
                            st_rel.append(site.species_string, site.frac_coords)
                            
                    st_rel.to("POSCAR", "structure.{}.vasp".format(i))
                    energy_corrlist.append([energy_ref[i], energy])

                    # Examine relaxation
                    dmax = 0
                    dmax_id = 0
                    for j, site in enumerate(stbak):
                        d = site.distance(st_rel[j])
                        if d > dmax:
                            dmax = d
                            dmax_id = j
                    relax_max.append([dmax_id, dmax])

                else:
                    solver_input.update_info_by_structure(st_rel)
                    solver_input.write_input(os.path.join(fmtstr.format(i),"baseinput{}".format(runstep+1)))
                    rundir_list.append(
                        os.path.abspath(
                            os.path.join(fmtstr.format(i),"baseinput{}".format(runstep+1)
                            )
                        )
                    )
            rundir_list = comm.gather(rundir_list, root = 0)
            if myreplica == 0 and not finalrun:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                with open(os.path.join(rootdir, "rundirs.txt"), "w") as fi:
                    fi.write("\n".join(rundir_list))
                    fi.flush()
                    os.fsync(fi.fileno())

                        
            if finalrun:
                np.savetxt("energy_corr.dat", energy_corrlist)
                with open("relax_max.dat", "w") as fi:
                    for row in relax_max:
                        fi.write("{} \t {}\n".format(row[0], row[1]))
                os.chdir(rootdir)
                comm.Barrier()
                if myreplica == 0:
                    with open("ALloop.progress", "a") as fi:
                        fi.write("AL{}\n".format(ALstep))
                        fi.flush()
                        os.fsync(fi.fileno())
                sys.exit(0)
            else:
                comm.Barrier()
                if myreplica == 0:
                    with open(os.path.join(os.path.pardir, "baseinput.progress"), "a") as fi:
                        fi.write(str(runstep) + "\n")
                        fi.flush()
                        os.fsync(fi.fileno())
                sys.exit(0)


        else: # Create first input for this AL step
            os.makedirs(ALdir, exist_ok=False)
            solver_input.from_directory(alparams.base_input_dir[0])
            os.chdir(ALdir)
            rundir_list = []
            for i in range(0, nsteps, sample_frequency):
                # In st_in, used for aenet,
                # (i)  "Ignore species" are omitted and should be added.
                # (ii) "Vacancy element" is included and should be removed.

                st_in = Structure.from_file(
                    os.path.join(MCdir, str(myreplica), "structure.{}.vasp".format(i))
                )

                # (i)
                st = map2perflat(perf_st, st_in)
                st.remove_species(["X"])
                # (ii)
                if alparams.vac_space_holder:
                    st.remove_species(alparams.vac_space_holder)
                perturb_structure(st,perturb)
                solver_input.update_info_by_structure(st)
                solver_input.write_input(
                    os.path.join(fmtstr.format(i), "baseinput0")
                )
                rundir_list.append(
                    os.path.abspath(
                        os.path.join(
                            fmtstr.format(i), "baseinput0"
                        )
                    )
                )
            rundir_list = comm.gather(rundir_list, root = 0)
            if myreplica == 0:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                with open(os.path.join(rootdir, "rundirs.txt"), "w") as fi:
                    fi.write("\n".join(rundir_list))
                    fi.flush()
                    os.fsync(fi.fileno())

def main():
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile)


if __name__ == "__main__":
    main()
