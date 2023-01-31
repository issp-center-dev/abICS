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

import logging
logger = logging.getLogger("main")


def main_potts(params_root: MutableMapping):
    param_config = params_root["config"]
    Q = param_config.get("Q", 2)
    Ls = param_config["L"]
    nspins = typing.cast(int, np.product(Ls))
    write_node = True

    model = Potts()
    sampler_type = params_root["sampling"].get("sampler", "RXMC")
    if sampler_type == "RXMC":
        rxparams = RXParams.from_dict(params_root["sampling"])
        nreplicas = rxparams.nreplicas
        configs = [Configuration(Q, Ls) for _ in range(nreplicas)]

        comm = RX_MPI_init(rxparams.nreplicas, rxparams.seed)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTs = np.linspace(kTstart, kTend, nreplicas)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        observer = Observer(comm, Lreload, {})

        nsteps = rxparams.nsteps
        RXtrial_frequency = rxparams.RXtrial_frequency
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency

        logger.info(f"-Running RXMC calculation with {nreplicas} replicas")
        logger.info(f"--Temperatures are linearly spaced from {kTstart} K to {kTend} K")

        RXcalc = TemperatureRX_MPI(
            comm, CanonicalMonteCarlo, model, configs, kTs, write_node=write_node
        )
        if Lreload:
            logger.info("-Reloading from previous calculation")
            RXcalc.reload()
        logger.info("-Starting RXMC calculation")
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
        pamcparams = PAMCParams.from_dict(params_root["sampling"])
        nreplicas = pamcparams.nreplicas
        configs = [Configuration(Q, Ls) for _ in range(nreplicas)]

        comm = RX_MPI_init(pamcparams.nreplicas, pamcparams.seed)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = pamcparams.kTstart
        kTend = pamcparams.kTend
        kTnum = pamcparams.kTnum
        kTs = np.linspace(kTstart, kTend, kTnum)

        # Set Lreload to True when restarting
        Lreload = pamcparams.reload

        observer = Observer(comm, Lreload, {})

        nsteps = pamcparams.nsteps
        resample_frequency = pamcparams.resample_frequency
        sample_frequency = pamcparams.sample_frequency
        print_frequency = pamcparams.print_frequency

        logger.info(f"-Running PAMC calculation with {nreplicas} replicas")
        logger.info(f"--Temperatures are linearly spaced from {kTstart} K to {kTend} K")

        calc = PopulationAnnealing(
            comm, CanonicalMonteCarlo, model, configs, kTs, write_node=write_node
        )
        if Lreload:
            logger.info("-Reloading from previous calculation")
            calc.reload()
        logger.info("-Starting PAMC calculation")
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
        comm = RX_MPI_init(rxparams.nreplicas, rxparams.seed)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        observer = Observer(comm, Lreload, {})

        nsteps = rxparams.nsteps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency
        logger.info(f"-Running parallel random sampling")
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

        comm = RX_MPI_init(rxparams.nreplicas, rxparams.seed)

        # RXMC parameters
        # specify temperatures for each replica, number of steps, etc.
        kTstart = rxparams.kTstart
        kTend = rxparams.kTend
        kTs = np.linspace(kTstart, kTend, nreplicas)

        # Set Lreload to True when restarting
        Lreload = rxparams.reload

        observer = Observer(comm, Lreload, {})

        nsteps = rxparams.nsteps
        sample_frequency = rxparams.sample_frequency
        print_frequency = rxparams.print_frequency
        logger.info(f"-Running parallel MC sampling")
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
        logger.error("Unknown sampler. Exiting...")
        sys.exit(1)

    logger.info("--Sampling completed sucessfully.")
    logger.info("Exiting normally on {}\n".format(datetime.datetime.now()))
