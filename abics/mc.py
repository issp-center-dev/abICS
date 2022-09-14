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

from typing import Any, List, Tuple, TextIO, Union
from numpy.typing import NDArray

from math import exp

from abc import ABCMeta, abstractmethod

# import SFMT_cython.sfmt_random as sfmt_random
# from multiprocessing import Process, Queue, Pool, TimeoutError
import os
import sys
import numpy as np
import numpy.random as rand

from abics import __version__

verylargeint = sys.maxsize


"""Defines base classes for Monte Carlo simulations"""


class Model(metaclass=ABCMeta):
    """This class defines a model whose energy equals 0 no matter the configuration, and the configuration
    never changes.
    This is a base template for building useful models."""

    model_name = ""

    # def __init__(self):

    @abstractmethod
    def energy(self, config) -> float:
        """
        Calculate energy of configuration: input: config

        Parameters
        ----------
        config: config object
            configuration

        Returns
        -------
        energy: float

        """
        ...

    @abstractmethod
    def trialstep(self, config, energy: float) -> Tuple[Any, float]:
        """Define a trial step on config

        Returns dconfig, which can contain the minimal information for
        constructing the trial configuration from config to be used in newconfig().
        Make sure that config is unchanged.

        Parameters
        ----------
        config: config object
            current configuration

        energy: float
            current energy

        Returns
        -------
        dconfig: config object
            The minimal information for constructing the trial configuration
            from config to be used in newconfig()
        dE: float
            Energy difference
        """

        # Return only change in configuration dconfig so that
        # you don't have to copy entire configurations,
        # which can sometimes be costly
        ...

    @abstractmethod
    def newconfig(self, config, dconfig):
        """
        Update config by using the trial step, dconfig

        Parameters
        ----------
        config: config object
            Original configuration
            This may be mutated through this function
        dconfig: config object
            Difference of configuration

        Returns
        -------
        config: config object
            updated configuration
        """
        return config


class Grid1D:
    def __init__(self, dx, minx, maxx):
        """

        Parameters
        ----------
        dx: float
            interval
        minx: float
            minimum value of x
        maxx: float
            maximum value of x

        """
        self.dx = dx
        self.x = np.arange(minx, maxx, dx)


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


class ObsInfo:
    nargs: int
    lengths: List[int]

    def __init__(self, *args):
        """

        Parameters
        ----------
        args: list
        """
        self.nargs = len(args)
        self.lengths = []
        for arg in args:
            # Inelegant way to make everything a 1D array
            arg = np.array([arg])
            arg = arg.ravel()
            self.lengths.append(len(arg))

    def decode(self, obs_array):
        """

        Parameters
        ----------
        obs_array: numpy array

        Returns
        -------
        obs: list

        """

        obs = []
        idx = 0
        for i in range(self.nargs):
            length = self.lengths[i]
            if length == 1:
                obs.append(obs_array[idx])
            else:
                obs.append(obs_array[idx : idx + length])
            idx += length
        return obs


# @profile
def obs_encode(*args):
    """make pure 1D data

    Parameters
    ----------
    args: list

    Returns
    -------
    obs_array: numpy array

    """
    # nargs = np.array([len(args)])
    # args_length_list = []
    obs_array = np.array([])
    for arg in args:
        # Inelegant way to make everything a 1D array
        arg = np.array([arg])
        arg = arg.ravel()
        obs_array = np.concatenate((obs_array, arg))
        # args_length_list.append(len(arg))
    # args_length_array = np.array(args_length_list)
    # args_info = np.concatenate((nargs, args_length_array))
    return obs_array


class ObserverBase:
    def __init__(self):
        self.lprintcount = 0

    def obs_info(self, calc_state: "MCAlgorithm") -> ObsInfo:
        """

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm

        Returns
        -------
        args_info: numpy array

        """
        obs_log = self.logfunc(calc_state)
        if not isinstance(obs_log, tuple):
            obs_log = (obs_log,)
        obs_ND = []
        obs_save = self.savefunc(calc_state)
        if len(obs_save) > 0:
            if isinstance(obs_save, tuple):
                for obs in obs_save:
                    obs_ND.append(obs)
            else:
                obs_ND.append(obs_save)
        return ObsInfo(*obs_log, *obs_ND)

    def logfunc(self, calc_state: "MCAlgorithm") -> Tuple[float]:
        """returns values of observables

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm

        Returns
        -------
        calc_state.energy: tuple
            (Energy)
        """
        return (calc_state.energy,)

    def savefunc(self, calc_state: "MCAlgorithm") -> Union[Tuple[float], Tuple[()]]:
        """returns values of observables, which will be not printed in observe method.

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm

        Returns
        -------
        """
        return ()

    def writefile(self, calc_state: "MCAlgorithm") -> None:
        """

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm

        Returns
        -------

        """
        return None

    def observe(self, calc_state: "MCAlgorithm", outputfi: TextIO, lprint=True):
        """

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm
        outputfi: _io.TextIOWrapper
            TextIOWrapper for output
        lprint: boolean
            if true, log info is outputted to TextIOWrapper

        Returns
        -------
        obs_log: numpy array
            log information about observation
        """
        obs_log = np.atleast_1d(self.logfunc(calc_state))
        if lprint:
            line = f"{self.lprintcount}\t"
            for p in calc_state.parameters():
                line += f"{p}\t"
            for x in obs_log:
                line += f"{x}\t"
            line += "\n"
            outputfi.write(line)
            outputfi.flush()
            self.writefile(calc_state)
            self.lprintcount += 1
        obs_save = np.atleast_1d(self.savefunc(calc_state))
        if len(obs_save) > 0:
            obs_save = np.atleast_1d(obs_save)
            obs_save = obs_save.ravel()
            print(obs_log.shape, obs_save.shape)
            return np.concatenate((obs_log, obs_save))
        else:
            return obs_log


class MCAlgorithm(metaclass=ABCMeta):

    model: Model
    config: Any
    obs_save: List[NDArray]

    @abstractmethod
    def __init__(self, *args):
        ...

    @abstractmethod
    def MCstep(self) -> None:
        """perform one MC step"""
        ...

    @abstractmethod
    def parameters(self) -> List:
        """returns parameters (e.g., temperature)"""
        ...

    def run(
        self,
        nsteps: int,
        sample_frequency: int = verylargeint,
        print_frequency: int = verylargeint,
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
        observer: observer object
        save_obs: boolean

        Returns
        -------
        observables: list
        """

        observables = 0.0
        nsample = 0
        self.energy = self.model.energy(self.config)

        output = open("obs.dat", "a")
        with open(os.devnull, "w") as f:
            if hasattr(observer.observe(self, f), "__add__"):
                observe = True
            else:
                observe = False
        for i in range(1, nsteps + 1):
            self.MCstep()
            sys.stdout.flush()
            if observe and i % sample_frequency == 0:
                obs_step = observer.observe(self, output, i % print_frequency == 0)
                observables += obs_step
                nsample += 1
                if save_obs:
                    self.obs_save.append(obs_step)
        output.close()

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
        self.model = model
        self.config = config
        self.kT = kT
        self.obs_save = []

    # @profile
    def MCstep(self):
        dconfig, dE = self.model.trialstep(self.config, self.energy)
        # if self.energy == float("inf"):
        #    self.config = self.model.newconfig(self.config, dconfig)
        #    self.energy = dE
        accepted = True
        if dE >= 0.0:
            accept_probability = exp(-dE / self.kT)
            accepted = rand.random() <= accept_probability
        if accepted:
            self.config = self.model.newconfig(self.config, dconfig)
            self.energy += dE

    def parameters(self):
        return [self.kT]


class RandomSampling(CanonicalMonteCarlo):
    def MCstep(self):
        self.config.shuffle()
        self.energy = self.model.energy(self.config)


# For backward compatibility
if __version__ < "3":
    model = Model
    grid_1D = Grid1D
    observer_base = ObserverBase
