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

import copy
import sys,os

from mpi4py import MPI
import numpy as np
import scipy.constants as constants

from abics.mc import CanonicalMonteCarlo, RandomSampling
from abics.mc_mpi import (
    RX_MPI_init,
    TemperatureRX_MPI,
    RXParams,
    SamplerParams,
    ParallelRandomParams,
    EmbarrassinglyParallelSampling,
)
from abics.applications.latgas_abinitio_interface import (
    default_observer, 
    EnsembleParams,
    ensemble_error_observer,
)
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
    runner_ensemble,
    runner_multistep,
)
from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.qe import QESolver
from abics.applications.latgas_abinitio_interface.aenet import aenetSolver
from abics.applications.latgas_abinitio_interface.openmx import OpenMXSolver
from abics.applications.latgas_abinitio_interface.mocksolver import MockSolver
from abics.applications.latgas_abinitio_interface.params import DFTParams

from abics.util import exists_on_all_nodes

def main_impl(tomlfile, ALrun=False):
    samplerparams = SamplerParams.from_toml(tomlfile)
    if samplerparams.sampler == "RXMC":
        rxparams = RXParams.from_toml(tomlfile)
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica

        kB = constants.value(u"Boltzmann constant in eV/K")

        comm = RX_MPI_init(rxparams)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTs = kB * np.linspace(kTstart, kTend, nreplicas)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        RXtrial_frequency = rxparams.RXtrial_frequency
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency

    elif samplerparams.sampler == "parallelRand":
        rxparams = ParallelRandomParams.from_toml(tomlfile)
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica
        comm = RX_MPI_init(rxparams)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency

    elif samplerparams.sampler == "parallelMC":
        rxparams = RXParams.from_toml(tomlfile)
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica

        kB = constants.value(u"Boltzmann constant in eV/K")

        comm = RX_MPI_init(rxparams)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTs = kB * np.linspace(kTstart, kTend, nreplicas)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency

    else:
        print("Unknown sampler. Exiting...")
        sys.exit(1)

    dftparams = DFTParams.from_toml(tomlfile)

    if dftparams.solver == "vasp":
        solver = VASPSolver(dftparams.path)
    elif dftparams.solver == "qe":
        parallel_level = dftparams.properties.get("parallel_level", {})
        solver = QESolver(dftparams.path, parallel_level=parallel_level)
    elif dftparams.solver == "aenet":
        solver = aenetSolver(
            dftparams.path, dftparams.ignore_species, dftparams.solver_run_scheme
        )
    elif dftparams.solver == "openmx":
        solver = OpenMXSolver(dftparams.path)
    elif dftparams.solver == "mock":
        solver = MockSolver()
    else:
        print("unknown solver: {}".format(dftparams.solver))
        sys.exit(1)

    # model setup
    # we first choose a "model" defining how to perform energy calculations and trial steps
    # on the "configuration" defined below
    if dftparams.ensemble:
        if len(dftparams.base_input_dir) ==1:
            print("You must specify more than one base_input_dir for ensemble calculator")
            sys.exit(1)
        energy_calculator = runner_ensemble(
            base_input_dirs=dftparams.base_input_dir,
            Solver=solver,
            runner=runner,
            nprocs_per_solver=nprocs_per_replica,
            comm=MPI.COMM_SELF,
            perturb=dftparams.perturb,
            solver_run_scheme=dftparams.solver_run_scheme,
            use_tmpdir=dftparams.use_tmpdir,
        )
    else:
        if len(dftparams.base_input_dir) == 1:
            energy_calculator = runner(
                base_input_dir=dftparams.base_input_dir[0],
                Solver=solver,
                nprocs_per_solver=nprocs_per_replica,
                comm=MPI.COMM_SELF,
                perturb=dftparams.perturb,
                solver_run_scheme=dftparams.solver_run_scheme,
                use_tmpdir=dftparams.use_tmpdir,
            )
        else:
            energy_calculator = runner_multistep(
                base_input_dirs=dftparams.base_input_dir,
                Solver=solver,
                runner=runner,
                nprocs_per_solver=nprocs_per_replica,
                comm=MPI.COMM_SELF,
                perturb=dftparams.perturb,
                solver_run_scheme=dftparams.solver_run_scheme,
                use_tmpdir=dftparams.use_tmpdir,
            )
    model = dft_latgas(energy_calculator, save_history=False)

    # defect sublattice setup

    configparams = DFTConfigParams.from_toml(tomlfile)

    spinel_config = defect_config(configparams)

    # configs = []
    # for i in range(nreplicas):
    #    configs.append(copy.deepcopy(spinel_config))
    configs = [spinel_config] * nreplicas

    obsparams = ObserverParams.from_toml(tomlfile)

    # NNP ensemble error estimation
    ensembleparams = EnsembleParams.from_toml(tomlfile)
    if ensembleparams:
        if ensembleparams.solver == "vasp":
            solver = VASPSolver(ensembleparams.path)
        elif ensembleparams.solver == "qe":
            parallel_level = ensembleparams.properties.get("parallel_level", {})
            solver = QESolver(ensembleparams.path, parallel_level=parallel_level)
        elif ensembleparams.solver == "aenet":
            solver = aenetSolver(
                ensembleparams.path, ensembleparams.ignore_species, ensembleparams.solver_run_scheme
            )
        elif ensembleparams.solver == "openmx":
            solver = OpenMXSolver(ensembleparams.path)
        elif ensembleparams.solver == "mock":
            solver = MockSolver()
        else:
            print("unknown solver: {}".format(ensembleparams.solver))
            sys.exit(1)
    
        energy_calculators = [
            runner(
                base_input_dir=base_input_dir,
                Solver=copy.deepcopy(solver),
                nprocs_per_solver=nprocs_per_replica,
                comm=MPI.COMM_SELF,
                perturb=ensembleparams.perturb,
                solver_run_scheme=ensembleparams.solver_run_scheme,
            ) for base_input_dir in ensembleparams.base_input_dirs
        ]
        print(ensembleparams.base_input_dirs)
        observer = ensemble_error_observer(comm, energy_calculators, Lreload)
    else:
        observer = default_observer(comm, Lreload)

    # Active learning mode
    if ALrun:
        if "train" in os.listdir():
            # Check how many AL iterations have been performed
            i = 0
            while exists_on_all_nodes(comm, "MC{}".format(i)):
                i += 1
            with open("ALloop.progress", "r") as fi:
                last_li = fi.readlines(-1)[-1]
            if "train" not in last_li:
                print("You should train before next MC sampling.")
                sys.exit(1)
            if Lreload:
                rootdir = os.getcwd()
                os.chdir("MC{}".format(i-1))
                MCid = i - 1
            else:
                # Make new directory and perform sampling there
                if comm.Get_rank() == 0:
                    os.mkdir("MC{}".format(i))
                comm.Barrier()
                rootdir = os.getcwd()
                os.chdir("MC{}".format(i))
                MCid = i
        else:
            print("You should train before MC sampling in AL mode.")
            sys.exit(1)


    if samplerparams.sampler == "RXMC":
        # RXMC calculation
        RXcalc = TemperatureRX_MPI(comm, CanonicalMonteCarlo, model, configs, kTs)
        if Lreload:
            RXcalc.reload()
        obs = RXcalc.run(
            nsteps,
            RXtrial_frequency,
            sample_frequency,
            print_frequency,
            observer=observer,
            subdirs=True,
        )

        if comm.Get_rank() == 0:
            print(obs)

    elif samplerparams.sampler == "parallelRand":
        calc = EmbarrassinglyParallelSampling(comm, RandomSampling, model, configs)
        if Lreload:
            calc.reload()
        obs = calc.run(
            nsteps,
            sample_frequency,
            print_frequency,
            observer=observer,
            subdirs=True,
        )

        if comm.Get_rank() == 0:
            print(obs)

    elif samplerparams.sampler == "parallelMC":
        calc = EmbarrassinglyParallelSampling(
            comm, CanonicalMonteCarlo, model, configs, kTs
        )
        if Lreload:
            calc.reload()
        obs = calc.run(
            nsteps,
            sample_frequency,
            print_frequency,
            observer=observer,
            subdirs=True,
        )

        if comm.Get_rank() == 0:
            print(obs)
    if ALrun:
        os.chdir(rootdir)
        if comm.Get_rank() == 0:
            with open("ALloop.progress", "a") as fi:
                fi.write("MC{}\n".format(MCid))
                fi.flush()
                os.fsync(fi.fileno())

def main():
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile)

def mainAL():
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile, True)



if __name__ == "__main__":
    main()
