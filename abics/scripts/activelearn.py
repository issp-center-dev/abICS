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
import os

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
        solver = VASPSolver(alparams.path)
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

    # model setup
    # we first choose a "model" defining how to perform energy calculations and trial steps
    # on the "configuration" defined below
    if len(alparams.base_input_dir) == 1:
        energy_calculator = runner(
            base_input_dir=alparams.base_input_dir[0],
            Solver=solver,
            nprocs_per_solver=nprocs_per_replica,
            comm=MPI.COMM_SELF,
            perturb=alparams.perturb,
            solver_run_scheme=alparams.solver_run_scheme,
        )
    else:
        energy_calculator = runner_multistep(
            base_input_dirs=alparams.base_input_dir,
            Solver=solver,
            runner=runner,
            nprocs_per_solver=nprocs_per_replica,
            comm=MPI.COMM_SELF,
            perturb=alparams.perturb,
            solver_run_scheme=alparams.solver_run_scheme,
        )

    myreplica = comm.Get_rank()

    # Find newest MC run
    i = 0
    while os.path.exists("MC{}".format(i)):
        i += 1

    rootdir = os.getcwd()
    if i == 0:  # Random sampling!
        if os.path.exists("AL0"):
            print(
                "It seems you've already run the first active learning step. You should train now."
            )
            sys.exit(1)
        model = dft_latgas(energy_calculator, save_history=False)
        configparams = DFTConfigParams.from_toml(tomlfile)
        config = defect_config(configparams)
        configs = [config] * nreplicas
        if myreplica == 0:
            os.mkdir("AL0")
        comm.Barrier()

        os.chdir("AL0")

        calc = EmbarrassinglyParallelSampling(comm, RandomSampling, model, configs)
        obs = calc.run(
            nsteps=nsteps // sample_frequency,
            sample_frequency=1,
            print_frequency=1,
            observer=default_observer(comm, False),
            subdirs=True,
        )

        if comm.Get_rank() == 0:
            print(obs)
        os.chdir(rootdir)

        if comm.Get_rank() == 0:
            with open("ALloop.progress", "a") as fi:
                fi.write("AL0\n")

    else:  # Active learning!
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
        os.makedirs(ALdir, exist_ok=False)
        os.chdir(ALdir)
        energy_corrlist = []
        relax_max = []
        config = defect_config(configparams)
        perf_st = config.structure
        if alparams.ignore_species:
            ignore_structure = perf_st.copy()
            remove_sp = filter(
                lambda sp: sp not in alparams.ignore_species, perf_st.symbol_set
            )
            ignore_structure.remove_species(remove_sp)

        for i in range(0, nsteps, sample_frequency):
            # In st_in, used for aenet,
            # (i)  "Ignore species" are omitted and should be added.
            # (ii) "Vacancy element" is included and should be removed.
            st_in = Structure.from_file(
                os.path.join(MCdir, str(myreplica), "structure.{}.vasp".format(i))
            )

            # (i)
            st = map2perflat(perf_st, st_in)
            # (ii)
            if alparams.vac_space_holder:
                st.remove_species(alparams.vac_space_holder)
            stbak = st.copy()

            energy, st_rel = energy_calculator.submit(st, os.path.join(ALdir, "output"))

            # energy_calculator may return the structure w/o ignored structure
            if alparams.ignore_species:
                for site in ignore_structure:
                    st_rel.append(site.species_string, site.frac_coords)
            st_rel.to("POSCAR", "structure.{}.vasp".format(i))
            energy_corrlist.append([energy_ref[i], energy])
            np.savetxt("energy_corr.dat", energy_corrlist)

            dmax = 0
            dmax_id = 0
            for j, site in enumerate(stbak):
                d = site.distance(st_rel[j])
                if d > dmax:
                    dmax = d
                    dmax_id = j
            relax_max.append([dmax_id, dmax])
            with open("relax_max.dat", "w") as fi:
                for row in relax_max:
                    fi.write("{} \t {}\n".format(row[0], row[1]))
        os.chdir(rootdir)
        if comm.Get_rank() == 0:
            with open("ALloop.progress", "a") as fi:
                fi.write("AL{}\n".format(ALstep))


def main():
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile)


if __name__ == "__main__":
    main()
