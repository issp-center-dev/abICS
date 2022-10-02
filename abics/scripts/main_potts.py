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

import typing
from typing import MutableMapping

import sys
import datetime

import numpy as np

from abics import __version__
from abics.mc import CanonicalMonteCarlo, RandomSampling
from abics.sampling.mc_mpi import RX_MPI_init
from abics.sampling.rxmc import TemperatureRX_MPI, RXParams
from abics.sampling.pamc import PopulationAnnealing, PAMCParams
from abics.sampling.simple_parallel import EmbarrassinglyParallelSampling, ParallelRandomParams

from abics.applications.lattice_model.potts import Potts, Configuration, Observer


def main_potts(params_root: MutableMapping):
    param_config = params_root["config"]
    Q = param_config.get("Q", 2)
    Ls = param_config["L"]
    nspins = typing.cast(int, np.product(Ls))
    write_node = True

    model = Potts()
    observer = Observer()
    sampler_type = params_root["sampling"].get("sampler", "RXMC")
    if sampler_type == "RXMC":
        rxparams = RXParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        configs = [Configuration(Q, Ls) for _ in range(nreplicas)]

        comm = RX_MPI_init(rxparams)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTs = np.linspace(kTstart, kTend, nreplicas)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        RXtrial_frequency = rxparams.RXtrial_frequency
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency

        if comm.Get_rank() == 0:
            print(f"-Running RXMC calculation with {nreplicas} replicas")
            print(f"--Temperatures are linearly spaced from {kTstart} K to {kTend} K")
            sys.stdout.flush()

        RXcalc = TemperatureRX_MPI(
            comm, CanonicalMonteCarlo, model, configs, kTs, write_node=write_node
        )
        if Lreload:
            if comm.Get_rank() == 0:
                print("-Reloading from previous calculation")
            RXcalc.reload()
        if comm.Get_rank() == 0:
            print("-Starting RXMC calculation")
            sys.stdout.flush()
        obs = RXcalc.run(
            nsteps,
            RXtrial_frequency=RXtrial_frequency,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            nsubsteps_in_step=nspins,
            observer=observer,
            subdirs=True,
        )

    elif sampler_type == "PAMC":
        rxparams = PAMCParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        configs = [Configuration(Q, Ls) for _ in range(nreplicas)]

        comm = RX_MPI_init(rxparams)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTnum = rxparams.kTnum
        kTs = np.linspace(kTstart, kTend, kTnum)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        resample_frequency = rxparams.resample_frequency
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency

        if comm.Get_rank() == 0:
            print(f"-Running RXMC calculation with {nreplicas} replicas")
            print(f"--Temperatures are linearly spaced from {kTstart} K to {kTend} K")
            sys.stdout.flush()

        calc = PopulationAnnealing(
            comm, CanonicalMonteCarlo, model, configs, kTs, write_node=write_node
        )
        if Lreload:
            if comm.Get_rank() == 0:
                print("-Reloading from previous calculation")
            calc.reload()
        if comm.Get_rank() == 0:
            print("-Starting RXMC calculation")
            sys.stdout.flush()
        obs = calc.run(
            nsteps,
            resample_frequency=resample_frequency,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            nsubsteps_in_step=nspins,
            observer=observer,
            subdirs=True,
        )

    elif sampler_type == "parallelRand":
        rxparams = ParallelRandomParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        comm = RX_MPI_init(rxparams)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency
        if comm.Get_rank() == 0:
            print(f"-Running parallel random sampling")
            sys.stdout.flush()
        calc = EmbarrassinglyParallelSampling(
            comm, RandomSampling, model, configs, write_node=write_node
        )
        if Lreload:
            calc.reload()
        obs = calc.run(
            nsteps,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            nsubsteps_in_step=nspins,
            observer=observer,
            subdirs=True,
        )
    elif sampler_type == "parallelMC":
        rxparams = RXParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas

        comm = RX_MPI_init(rxparams)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTs = np.linspace(kTstart, kTend, nreplicas)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        nsteps = rxparams.nsteps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency
        if comm.Get_rank() == 0:
            print(f"-Running parallel MC sampling")
            sys.stdout.flush()
        calc = EmbarrassinglyParallelSampling(
            comm, CanonicalMonteCarlo, model, configs, kTs, write_node=write_node
        )
        if Lreload:
            calc.reload()
        obs = calc.run(
            nsteps,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            nsubsteps_in_step=nspins,
            observer=observer,
            subdirs=True,
        )
    else:
        print("Unknown sampler. Exiting...")
        sys.exit(1)

    if comm.Get_rank() == 0:
        print("--Sampling completed sucessfully.")

    if comm.Get_rank() == 0:
        now = datetime.datetime.now()
        print(f"Exiting normally on {now}\n")
