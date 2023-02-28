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

from __future__ import annotations

from typing import Any, Optional

from abc import ABCMeta, abstractmethod

# import SFMT_cython.sfmt_random as sfmt_random
# from multiprocessing import Process, Queue, Pool, TimeoutError
import os
import sys
import numpy as np
import numpy.random as rand

from abics import __version__
from abics.observer import ObserverBase
from abics.model import Model

import logging
logger = logging.getLogger("main")

verylargeint = sys.maxsize

def binning(x, nlevels: int):
    """

    Parameters
    ----------
    x: list
        Coordinates
    nlevels: int
        Number to determine the number of measurements

    Returns
    -------
    error_estimate: list
        error estimation
    """
    error_estimate = []
    x = np.array(x, dtype=np.float64)
    # assert 2**nlevels*10 < len(x)
    throwout = len(x) % (2**nlevels)
    if throwout != 0:
        # The number of measurements must be divisible by 2**nlevels
        # If not, throw out initial measurements
        x = x[throwout:]
    error_estimate.append(np.sqrt(np.var(x, ddof=1) / len(x)))
    for lvl in range(1, nlevels + 1):
        x_tmp = x
        x = (x_tmp[0::2] + x_tmp[1::2]) / 2.0
        error_estimate.append(np.sqrt(np.var(x, ddof=1) / len(x)))
    return error_estimate


class MCAlgorithm(metaclass=ABCMeta):

    model: Model
    config: Optional[Any]
    obs_save: list[np.ndarray]
    kT: float
    naccepted: int
    ntrials: int
    nfail: int

    @abstractmethod
    def __init__(self, *args):
        self.naccepted = 0
        self.ntrials = 0
        self.nfail = 0

    @abstractmethod
    def MCstep(self, nsubstep_in_step: int = 1) -> None:
        """perform one MC step"""
        ...

    @abstractmethod
    def parameters(self) -> list:
        """returns parameters (e.g., temperature)"""
        ...

    def run(
        self,
        nsteps: int,
        sample_frequency: int = verylargeint,
        print_frequency: int = verylargeint,
        nsubsteps_in_step: int = 1,
        observer: ObserverBase = ObserverBase(),
        save_obs: bool = False,
    ):
        """

        Parameters
        ----------
        nsteps: int
            The number of Monte Carlo steps for running.
        sample_frequency: int
            The number of Monte Carlo steps for observation of physical quantities.
        print_frequency: int
            The number of Monte Carlo steps for saving physical quantities.
        nsubsteps_in_step: int
            The number of Monte Carlo substeps in one MC step
        observer: observer object
        save_obs: boolean

        Returns
        -------
        observables: list
        """

        nsample = 0
        self.energy = self.model.energy(self.config)

        output = open("obs.dat", "a")
        nobs = 0
        with open(os.devnull, "w") as f:
            o = observer.observe(self, f)
            nobs = np.atleast_1d(o).flatten().size
        observables = np.zeros(nobs)
        for i in range(1, nsteps + 1):
            self.MCstep(nsubsteps_in_step)
            sys.stdout.flush()
            if i % sample_frequency == 0:
                obs_step = observer.observe(self, output, i % print_frequency == 0)
                observables += np.atleast_1d(obs_step).flatten()
                nsample += 1
                if save_obs:
                    self.obs_save.append(obs_step)
        output.close()

        logger.info("MCAlgorithm: acceptance = {}".format(self.naccepted/self.ntrials))
        logger.info("MCAlgotithm: nfail = {} / ntrials = {}".format(self.nfail, self.ntrials))

        if save_obs:
            np.save(open("obs_save.npy", "wb"), np.array(self.obs_save))
        if nsample > 0:
            observables /= nsample
            obs_info = observer.obs_info(self)
            return obs_info.decode(observables)
        else:
            return None


class CanonicalMonteCarlo(MCAlgorithm):
    def __init__(self, model: Model, kT: float, config):
        """

        Parameters
        ----------
        model: Model
        kT: float
            Temperature
        config: config object
            Configuration
        """
        super().__init__()
        self.model = model
        self.config = config
        self.kT = kT
        self.obs_save = []

        logger.info("CanonicalMonteCarlo: parameters")
        logger.info("  kT = {}".format(self.kT))


    # @profile
    def MCstep(self, nsubsteps_in_step: int = 1):
        for istep in range(nsubsteps_in_step):
            dconfig, dE = self.model.trialstep(self.config, self.energy)
            # if self.energy == float("inf"):
            #    self.config = self.model.newconfig(self.config, dconfig)
            #    self.energy = dE
            accepted = True
            if dE >= 0.0:
                accept_probability = np.exp(-dE / self.kT)
                trial_rand = rand.random()
                accepted = trial_rand <= accept_probability

            if accepted:
                self.config = self.model.newconfig(self.config, dconfig)
                self.energy += dE
                self.naccepted += 1
            self.ntrials += 1

            logger.debug("MCstep: dE/kT = {:e}".format(dE/self.kT))
            logger.debug("MCstep: dG    = {:e}".format(dG))
            if dG > 0.0:
                logger.debug("MCstep: prob  = {}".format(accept_probability))
                logger.debug("MCstep: rand  = {}".format(trial_rand))
            logger.debug("MCstep: {}".format("accepted" if accepted else "rejected"))

    def parameters(self):
        return [self.kT]

class WeightedCanonicalMonteCarlo(CanonicalMonteCarlo):
    def __init__(self, model: Model, kT: float, config):
        super().__init__(model, kT, config)

    # @profile
    def MCstep(self, nsubsteps_in_step: int = 1):
        for istep in range(nsubsteps_in_step):
            dconfig, dE, *opt = self.model.trialstep(self.config, self.energy)

            # ln weight factor
            dW = opt[0] if len(opt) > 0 else 0.0

            if np.isnan(dW):
                dG = float("nan")
                accepted = False
                self.nfail += 1
            else:
                dG = dE/self.kT - dW

                if dG >= 0.0:
                    accept_probability = np.exp(-dG)
                    trial_rand = rand.random()
                    accepted = trial_rand <= accept_probability
                else:
                    accepted = True

            if accepted:
                self.config = self.model.newconfig(self.config, dconfig)
                self.energy += dE
                self.naccepted += 1

            self.ntrials += 1

            logger.debug("MCstep: dE/kT = {:e}".format(dE/self.kT))
            logger.debug("MCstep: dG    = {:e}".format(dG))
            if dG >= 0.0:
                logger.debug("MCstep: prob  = {}".format(accept_probability))
                logger.debug("MCstep: rand  = {}".format(trial_rand))
            logger.debug("MCstep: {}".format("accepted" if accepted else "rejected"))

class RandomSampling(CanonicalMonteCarlo):
    def MCstep(self, nsubsteps_in_step: int = 1):
        if self.config is None:
            raise RuntimeError("config is not set")
        self.config.shuffle()
        self.energy = self.model.energy(self.config)
        self.ntrials += 1
        self.naccepted += 1
