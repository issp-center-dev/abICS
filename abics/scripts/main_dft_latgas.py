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

from typing import MutableMapping, Union

from mpi4py import MPI
import copy
import sys, os, shutil
import datetime

import numpy as np
import scipy.constants as constants

from abics import __version__
from abics.mc import CanonicalMonteCarlo, RandomSampling

from abics.sampling.mc_mpi import RX_MPI_init
from abics.sampling.rxmc import TemperatureRX_MPI, RXParams
from abics.sampling.pamc import PopulationAnnealing, PAMCParams
from abics.sampling.simple_parallel import EmbarrassinglyParallelSampling, ParallelRandomParams

from abics.applications.latgas_abinitio_interface import (
    DefaultObserver,
    EnsembleParams,
    EnsembleErrorObserver,
)
from abics.applications.latgas_abinitio_interface.model_setup import (
    DFTLatticeGas,
    ObserverParams,
)
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)
from abics.applications.latgas_abinitio_interface.run_base_mpi import (
    Runner,
    RunnerEnsemble,
    RunnerMultistep,
)
from abics.applications.latgas_abinitio_interface.base_solver import SolverBase
from abics.applications.latgas_abinitio_interface.params import DFTParams

from abics.util import exists_on_all_nodes

import logging
logger = logging.getLogger("main")

def main_dft_latgas(params_root: MutableMapping):
    dftparams = DFTParams.from_dict(params_root["sampling"]["solver"])
    sampler_type = params_root["sampling"].get("sampler", "RXMC")
    if sampler_type == "RXMC":
        rxparams = RXParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica

        kB = constants.value("Boltzmann constant in eV/K")

        nensemble = len(dftparams.base_input_dir)
        comm, commEnsemble, commAll = RX_MPI_init(rxparams.nreplicas, rxparams.seed, nensemble)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTs = kB * rxparams.kTs
        kTstart = rxparams.kTs[0]
        kTend = rxparams.kTs[-1]

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        RXtrial_frequency = rxparams.RXtrial_frequency
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency

        logger.info(f"-Running RXMC calculation with {nreplicas} replicas")
        logger.info(f"--Temperature varies from {kTstart} K to {kTend} K")

    elif sampler_type == "PAMC":
        pamcparams = PAMCParams.from_dict(params_root["sampling"])
        nreplicas = pamcparams.nreplicas
        nprocs_per_replica = pamcparams.nprocs_per_replica

        kB = constants.value("Boltzmann constant in eV/K")

        nensemble = len(dftparams.base_input_dir)
        comm, commEnsemble, commAll = RX_MPI_init(pamcparams.nreplicas, pamcparams.seed, nensemble)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = pamcparams.kTs[0]
        kTend = pamcparams.kTs[-1]
        if kTstart < kTend:
            kTstart, kTend = kTend, kTstart
        kTs = kB * pamcparams.kTs

        # Set Lreload to True when restarting
        Lreload = pamcparams.reload

        nsteps = pamcparams.nsteps
        resample_frequency = pamcparams.resample_frequency
        sample_frequency = pamcparams.sample_frequency
        print_frequency = pamcparams.print_frequency

        logger.info(f"-Running PAMC calculation with {nreplicas} replicas")
        logger.info(f"--Anneal from {kTstart} K to {kTend} K")

    elif sampler_type == "parallelRand":
        rxparams = ParallelRandomParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica
        nensemble = len(dftparams.base_input_dir)
        comm, commEnsemble, commAll = RX_MPI_init(rxparams.nreplicas, rxparams.seed, nensemble)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency
        logger.info(f"-Running parallel random sampling")

    elif sampler_type == "parallelMC":
        rxparams = RXParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica

        kB = constants.value("Boltzmann constant in eV/K")

        nensemble = len(dftparams.base_input_dir)
        comm, commEnsemble, commAll = RX_MPI_init(rxparams.nreplicas, rxparams.seed, nensemble)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTs[0]
        kTend = rxparams.kTs[-1]
        kTs = kB * rxparams.kTs

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency
        logger.info(f"-Running parallel MC sampling")

    else:
        logger.error("Unknown sampler. Exiting...")
        sys.exit(1)

    solvers = []
    for i in range(len(dftparams.base_input_dir)):
        solver: SolverBase = SolverBase.create(dftparams.solver, dftparams)
        solvers.append(solver)
    
    logger.info(f"-Setting up {dftparams.solver} solver for configuration energies")
    logger.info("--Base input is taken from {}".format(",".join(dftparams.base_input_dir)))

# >>>>>>> develop
    # model setup
    # we first choose a "model" defining how to perform energy calculations and trial steps
    # on the "configuration" defined below
    energy_calculator: Union[Runner, RunnerEnsemble, RunnerMultistep]
    if dftparams.ensemble:
        if len(dftparams.base_input_dir) == 1:
            logger.error("You must specify more than one base_input_dir for ensemble calculator")
            sys.exit(1)
        energy_calculator = RunnerEnsemble(
            base_input_dirs=dftparams.base_input_dir,
            Solvers=solvers,
            runner=Runner,
            nprocs_per_solver=nprocs_per_replica,
            comm=commEnsemble,
            perturb=dftparams.perturb,
            solver_run_scheme=dftparams.solver_run_scheme,
            use_tmpdir=dftparams.use_tmpdir,
        )
    else:
        if len(dftparams.base_input_dir) == 1:
            energy_calculator = Runner(
                base_input_dir=dftparams.base_input_dir[0],
                Solver=solvers[0],
                nprocs_per_solver=nprocs_per_replica,
                comm=MPI.COMM_SELF,
                perturb=dftparams.perturb,
                solver_run_scheme=dftparams.solver_run_scheme,
                use_tmpdir=dftparams.use_tmpdir,
            )
        else:
            energy_calculator = RunnerMultistep(
                base_input_dirs=dftparams.base_input_dir,
                Solvers=solvers,
                runner=Runner,
                nprocs_per_solver=nprocs_per_replica,
                comm=MPI.COMM_SELF,
                perturb=dftparams.perturb,
                solver_run_scheme=dftparams.solver_run_scheme,
                use_tmpdir=dftparams.use_tmpdir,
            )
    model = DFTLatticeGas(energy_calculator, save_history=False)
    logger.info("--Success.")
    logger.info("-Setting up the on-lattice model.")
    
    configparams = DFTConfigParams.from_dict(params_root["config"])

    spinel_config = defect_config(configparams)
    if configparams.init_structure is None:
        exit_code, msg = spinel_config.shuffle()
        logger.info("--Rank {}: {}".format(commAll.Get_rank(), msg))
        exit_code = commAll.allgather(exit_code)
        if np.sum(exit_code) != 0:
            logger.error("Failed to initialize configuration. Exiting.")
            sys.exit(1)
    else:
        logger.info("--Initial structure will be set to 'config.init_structure'.")
        spinel_config.reset_from_structure(configparams.init_structure)

    # configs = []
    # for i in range(nreplicas):
    #    configs.append(copy.deepcopy(spinel_config))
    configs = [spinel_config] * nreplicas

    obsparams = ObserverParams.from_dict(params_root["observer"])

    logger.info("--Success.")

    # NNP ensemble error estimation
    if "ensemble" in params_root:
        ensembleparams = EnsembleParams.from_dict(params_root["ensemble"])
        solver = SolverBase.create(ensembleparams.solver, ensembleparams)

        energy_calculators = [
            Runner(
                base_input_dir=base_input_dir,
                Solver=copy.deepcopy(solver),
                nprocs_per_solver=nprocs_per_replica,
                comm=MPI.COMM_SELF,
                perturb=ensembleparams.perturb,
                solver_run_scheme=ensembleparams.solver_run_scheme,
                use_tmpdir=dftparams.use_tmpdir,
            )
            for base_input_dir in ensembleparams.base_input_dirs
        ]
        observer: DefaultObserver = EnsembleErrorObserver(commEnsemble, energy_calculators, Lreload)
    else:
        observer = DefaultObserver(comm, Lreload)

    ALrun = exists_on_all_nodes(commAll, "ALloop.progress")

    # Active learning mode
    if ALrun:
        logger.info(f"-Running in active learning mode.")

        if "train0" in os.listdir():
            # Check how many AL iterations have been performed
            i = 0
            while exists_on_all_nodes(commAll, "MC{}".format(i)):
                i += 1
            with open("ALloop.progress", "r") as fi:
                last_li = fi.readlines(-1)[-1]
            if "train" not in last_li:
                logger.error("You should train before next MC sampling.")
                sys.exit(1)
            if Lreload:
                logger.info(f"--Restarting run in MC{i-1}")
                rootdir = os.getcwd()
                os.chdir("MC{}".format(i - 1))
                MCid = i - 1
            else:
                # Make new directory and perform sampling there
                if commAll.Get_rank() == 0:
                    logger.info(f"--MC sampling will be run in MC{i}")
                    os.mkdir("MC{}".format(i))
                    if dftparams.use_tmpdir:
                        logger.info(f"---Will use local tmpdir for {dftparams.solver} run")
                        # backup baseinput for this AL step
                        for j, d in enumerate(dftparams.base_input_dir):
                            shutil.copytree(d, "MC{}/baseinput{}".format(i, j))
                    sys.stdout.flush()
                commAll.Barrier()
                rootdir = os.getcwd()
                while not exists_on_all_nodes(commAll, "MC{}".format(i)):
                    pass
                os.chdir("MC{}".format(i))
                MCid = i
        else:
            logger.error("You should train before MC sampling in AL mode.")
            sys.exit(1)

    if commEnsemble.Get_rank() == 0:
        write_node = True
    else:
        write_node = False
    if sampler_type == "RXMC":
        # RXMC calculation
        RXcalc = TemperatureRX_MPI(
            comm, CanonicalMonteCarlo, model, configs, kTs, write_node=write_node
        )
        if Lreload:
            logger.info("-Reloading from previous calculation")
            RXcalc.reload()

        logger.info("-Starting RXMC calculation")
            
        obs = RXcalc.run(
            nsteps,
            RXtrial_frequency,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            observer=observer,
            subdirs=True,
        )

    elif sampler_type == "PAMC":
        # PAMC calculation
        PAcalc = PopulationAnnealing(
            comm, CanonicalMonteCarlo, model, configs, kTs, write_node=write_node
        )
        if Lreload:
            logger.info("-Reloading from previous calculation")
            PAcalc.reload()

        logger.info("-Starting PAMC calculation")
            
        obs = PAcalc.run(
            nsteps,
            resample_frequency,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            observer=observer,
            subdirs=True,
        )

    elif sampler_type == "parallelRand":
        calc = EmbarrassinglyParallelSampling(
            comm, RandomSampling, model, configs, write_node=write_node
        )
        if Lreload:
            calc.reload()
        obs = calc.run(
            nsteps,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            observer=observer,
            subdirs=True,
        )

    elif sampler_type == "parallelMC":
        calc = EmbarrassinglyParallelSampling(
            comm, CanonicalMonteCarlo, model, configs, kTs, write_node=write_node
        )
        if Lreload:
            calc.reload()
        obs = calc.run(
            nsteps,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            observer=observer,
            subdirs=True,
        )
    logger.info("--Sampling completed sucessfully.")

    if ALrun:
        os.chdir(rootdir)
        if comm.Get_rank() == 0 and write_node:
            logger.info("-Writing ALloop.progress")
            with open("ALloop.progress", "a") as fi:
                fi.write("MC{}\n".format(MCid))
                fi.flush()
                os.fsync(fi.fileno())
    logger.info("Exiting normally on {}\n".format(datetime.datetime.now()))
