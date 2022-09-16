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

from typing import MutableMapping

from mpi4py import MPI
import copy
import sys, os, shutil
import datetime

import numpy as np
import scipy.constants as constants

from abics import __version__
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
from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.qe import QESolver
from abics.applications.latgas_abinitio_interface.aenet import AenetSolver
from abics.applications.latgas_abinitio_interface.openmx import OpenMXSolver
from abics.applications.latgas_abinitio_interface.mocksolver import MockSolver
from abics.applications.latgas_abinitio_interface.params import DFTParams

from abics.util import exists_on_all_nodes


def main_dft_latgas(params_root: MutableMapping):
    dftparams = DFTParams.from_dict(params_root["sampling"]["solver"])
    sampler_type = params_root["sampling"].get("sampler", "RXMC")
    if sampler_type == "RXMC":
        rxparams = RXParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica

        kB = constants.value("Boltzmann constant in eV/K")

        nensemble = len(dftparams.base_input_dir)
        comm, commEnsemble, commAll = RX_MPI_init(rxparams, nensemble)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTs = kB * np.linspace(kTstart, kTend, nreplicas)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsweeps
        RXtrial_frequency = rxparams.RXtrial_frequency
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency

        if commAll.Get_rank() == 0:
            print(f"-Running RXMC calculation with {nreplicas} replicas")
            print(f"--Temperatures are linearly spaced from {kTstart} K to {kTend} K")
            sys.stdout.flush()

    elif sampler_type == "parallelRand":
        rxparams = ParallelRandomParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica
        nensemble = len(dftparams.base_input_dir)
        comm, commEnsemble, commAll = RX_MPI_init(rxparams, nensemble)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsweeps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency
        if commAll.Get_rank() == 0:
            print(f"-Running parallel random sampling")
            sys.stdout.flush()
    elif sampler_type == "parallelMC":
        rxparams = RXParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        nprocs_per_replica = rxparams.nprocs_per_replica

        kB = constants.value("Boltzmann constant in eV/K")

        nensemble = len(dftparams.base_input_dir)
        comm, commEnsemble, commAll = RX_MPI_init(rxparams, nensemble)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTs = kB * np.linspace(kTstart, kTend, nreplicas)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsweeps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency
        if commAll.Get_rank() == 0:
            print(f"-Running parallel MC sampling")
            sys.stdout.flush()
    else:
        print("Unknown sampler. Exiting...")
        sys.exit(1)

    if dftparams.solver == "vasp":
        solver = VASPSolver(dftparams.path)
    elif dftparams.solver == "qe":
        parallel_level = dftparams.properties.get("parallel_level", {})
        solver = QESolver(dftparams.path, parallel_level=parallel_level)
    elif dftparams.solver == "aenet":
        solver = AenetSolver(
            dftparams.path, dftparams.ignore_species, dftparams.solver_run_scheme
        )
    elif dftparams.solver == "openmx":
        solver = OpenMXSolver(dftparams.path)
    elif dftparams.solver == "mock":
        solver = MockSolver()
    else:
        print("unknown solver: {}".format(dftparams.solver))
        sys.exit(1)
    if commAll.Get_rank() == 0:
        print(f"-Setting up {dftparams.solver} solver for configuration energies")
        print("--Base input is taken from {}".format(",".join(dftparams.base_input_dir)))
        sys.stdout.flush()
    # model setup
    # we first choose a "model" defining how to perform energy calculations and trial steps
    # on the "configuration" defined below
    if dftparams.ensemble:
        if len(dftparams.base_input_dir) == 1:
            print(
                "You must specify more than one base_input_dir for ensemble calculator"
            )
            sys.exit(1)
        energy_calculator = RunnerEnsemble(
            base_input_dirs=dftparams.base_input_dir,
            Solver=solver,
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
                Solver=solver,
                nprocs_per_solver=nprocs_per_replica,
                comm=MPI.COMM_SELF,
                perturb=dftparams.perturb,
                solver_run_scheme=dftparams.solver_run_scheme,
                use_tmpdir=dftparams.use_tmpdir,
            )
        else:
            energy_calculator = RunnerMultistep(
                base_input_dirs=dftparams.base_input_dir,
                Solver=solver,
                runner=Runner,
                nprocs_per_solver=nprocs_per_replica,
                comm=MPI.COMM_SELF,
                perturb=dftparams.perturb,
                solver_run_scheme=dftparams.solver_run_scheme,
                use_tmpdir=dftparams.use_tmpdir,
            )
    model = DFTLatticeGas(energy_calculator, save_history=False)
    if commAll.Get_rank() == 0:
        print("--Success.")

    # defect sublattice setup
    if commAll.Get_rank() == 0:
        print("-Setting up the on-lattice model.")

    sys.stdout.flush()
    
    configparams = DFTConfigParams.from_dict(params_root["config"])

    spinel_config = defect_config(configparams)

    # configs = []
    # for i in range(nreplicas):
    #    configs.append(copy.deepcopy(spinel_config))
    configs = [spinel_config] * nreplicas

    obsparams = ObserverParams.from_dict(params_root["observer"])

    if commAll.Get_rank() == 0:
        print("--Success.")

    
    # NNP ensemble error estimation
    if "ensemble" in params_root:
        ensembleparams = EnsembleParams.from_dict(params_root["ensemble"])
        if ensembleparams.solver == "vasp":
            solver = VASPSolver(ensembleparams.path)
        elif ensembleparams.solver == "qe":
            parallel_level = ensembleparams.properties.get("parallel_level", {})
            solver = QESolver(ensembleparams.path, parallel_level=parallel_level)
        elif ensembleparams.solver == "aenet":
            solver = AenetSolver(
                ensembleparams.path,
                ensembleparams.ignore_species,
                ensembleparams.solver_run_scheme,
            )
        elif ensembleparams.solver == "openmx":
            solver = OpenMXSolver(ensembleparams.path)
        elif ensembleparams.solver == "mock":
            solver = MockSolver()
        else:
            print("unknown solver: {}".format(ensembleparams.solver))
            sys.exit(1)

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
        observer = EnsembleErrorObserver(commEnsemble, energy_calculators, Lreload)
    else:
        observer = DefaultObserver(comm, Lreload)

    ALrun = exists_on_all_nodes(commAll, "ALloop.progress")

    # Active learning mode
    if ALrun:
        if commAll.Get_rank() == 0:
            print(f"-Running in active learning mode.")

        if "train0" in os.listdir():
            # Check how many AL iterations have been performed
            i = 0
            while exists_on_all_nodes(commAll, "MC{}".format(i)):
                i += 1
            with open("ALloop.progress", "r") as fi:
                last_li = fi.readlines(-1)[-1]
            if "train" not in last_li:
                print("You should train before next MC sampling.")
                sys.exit(1)
            if Lreload:
                if commAll.Get_rank() == 0:
                    print(f"--Restarting run in MC{i-1}")
                    sys.stdout.flush()
                rootdir = os.getcwd()
                os.chdir("MC{}".format(i - 1))
                MCid = i - 1
            else:
                # Make new directory and perform sampling there
                if commAll.Get_rank() == 0:
                    print(f"--MC sampling will be run in MC{i}")
                    os.mkdir("MC{}".format(i))
                    if dftparams.use_tmpdir:
                        print(f"---Will use local tmpdir for {dftparams.solver} run")
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
            print("You should train before MC sampling in AL mode.")
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
            if commAll.Get_rank() == 0:
                print("-Reloading from previous calculation")
            RXcalc.reload()

        if commAll.Get_rank() == 0:
            print("-Starting RXMC calculation")
            sys.stdout.flush()
            
        obs = RXcalc.run(
            nsteps,
            RXtrial_frequency,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            observer=observer,
            subdirs=True,
        )

        #if comm.Get_rank() == 0 and write_node:
        #    print(obs)

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

        #if comm.Get_rank() == 0 and write_node:
        #    print(obs)

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

        #if comm.Get_rank() == 0 and write_node:
        #    print(obs)

    if commAll.Get_rank() == 0:
        print("--Sampling completed sucessfully.")
    if ALrun:
        os.chdir(rootdir)
        if comm.Get_rank() == 0 and write_node:
            print("-Writing ALloop.progress")
            with open("ALloop.progress", "a") as fi:
                fi.write("MC{}\n".format(MCid))
                fi.flush()
                os.fsync(fi.fileno())

    if commAll.Get_rank() == 0:
        now = datetime.datetime.now()
        print(f"Exiting normally on {now}\n")
