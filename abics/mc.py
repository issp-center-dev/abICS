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

from math import exp
from random import random

# import SFMT_cython.sfmt_random as sfmt_random
# from multiprocessing import Process, Queue, Pool, TimeoutError
import os
import sys
import numpy as np

verylargeint = sys.maxsize


"""Defines base classes for Monte Carlo simulations"""


class model:
    """ This class defines a model whose energy equals 0 no matter the configuration, and the configuration
    never changes.
    This is a base template for building useful models."""

    model_name = None

    # def __init__(self):

    def energy(self, config):
        """
        Calculate energy of configuration: input: config

        Parameters
        ----------
        config: config object
            configuration

        Returns
        -------

        """
        return 0.0

    def trialstep(self, config, energy):
        """
        Define a trial step on config. Returns dconfig, which can contain the minimal information for
        constructing the trial configuration from config to be used in newconfig().
        Make sure that config is unchanged.

        Parameters
        ----------
        config: config object
            configuration

        energy: float
            energy
        Returns
        -------
        dconfig: config object
            The minimal information for constructing the trial configuration from config to be used in newconfig()
        dE: float
            Energy diffence
        """
        dE = 0.0
        dconfig = None
        # Return only change in configuration dconfig so that
        # you don't have to copy entire configurations,
        # which can sometimes be costly
        return dconfig, dE

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


class grid_1D:
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


def binning(x, nlevels):
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
    throwout = len(x) % (2 ** nlevels)
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


empty_array = np.array([])
# @profile


def obs_encode(*args):
    """

    Parameters
    ----------
    args: list

    Returns
    -------
    obs_array: numpy array

    """
    # nargs = np.array([len(args)])
    # args_length_list = []
    obs_array = empty_array
    for arg in args:
        # Inelegant way to make everything a 1D array
        arg = np.array([arg])
        arg = arg.ravel()
        obs_array = np.concatenate((obs_array, arg))
        # args_length_list.append(len(arg))
    # args_length_array = np.array(args_length_list)
    # args_info = np.concatenate((nargs, args_length_array))
    return obs_array


def args_info(*args):
    """

    Parameters
    ----------
    args: list

    Returns
    -------
    args_info: numpy array
    """
    nargs = np.array([len(args)])
    args_length_list = []
    for arg in args:
        # Inelegant way to make everything a 1D array
        arg = np.array([arg])
        arg = arg.ravel()
        args_length_list.append(len(arg))
    args_length_array = np.array(args_length_list)
    args_info = np.concatenate((nargs, args_length_array))
    return args_info


def obs_decode(args_info, obs_array):
    """

    Parameters
    ----------
    args_info: numpy array
    obs_array: numpy array

    Returns
    -------
    args: list
    """
    nargs = args_info[0]
    args_length_array = args_info[1 : nargs + 1]
    args = []
    idx = 0
    for i in range(nargs):
        length = args_length_array[i]
        if length == 1:
            args.append(obs_array[idx])
        else:
            args.append(obs_array[idx : idx + length])
        idx += length
    return args


class observer_base:
    def __init__(self):
        self.lprintcount = 0

    def obs_info(self, calc_state):
        """

        Parameters
        ----------
        calc_state: MonteCarlo algorithm object
            MonteCarlo algorithm

        Returns
        -------
        args_info: numpy array

        """
        obs_log = self.logfunc(calc_state)
        if not isinstance(obs_log, tuple):
            obs_log = (obs_log,)
        obs_ND = []
        obs_save = self.savefuncs(calc_state)
        if obs_save is not None:
            if isinstance(obs_save, tuple):
                for obs in obs_save:
                    obs_ND.append(obs)
            else:
                obs_ND.append(obs_save)
        return args_info(*obs_log, *obs_ND)

    def logfunc(self, calc_state):
        """

        Parameters
        ----------
        calc_state: MonteCarlo algorithm object
            MonteCarlo algorithm

        Returns
        -------
        calc_state.energy: tuple
            (Energy)
        """
        return (calc_state.energy,)

    def savefuncs(self, calc_state):
        """

        Parameters
        ----------
        calc_state: MonteCarlo algorithm object
            MonteCarlo algorithm

        Returns
        -------
        """
        return None

    def writefile(self, calc_state):
        """

        Parameters
        ----------
        calc_state: MonteCarlo algorithm object
            MonteCarlo algorithm

        Returns
        -------

        """
        return None

    def observe(self, calc_state, outputfi, lprint=True):
        """

        Parameters
        ----------
        calc_state: MonteCarlo algorithm object
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
            outputfi.write(
                str(self.lprintcount)
                + "\t"
                + str(calc_state.kT)
                + "\t"
                + "\t".join([str(x) for x in obs_log])
                + "\n"
            )
            outputfi.flush()
            self.writefile(calc_state)
            self.lprintcount += 1
        obs_save = self.savefuncs(calc_state)
        if obs_save is not None:
            obs_save = np.atleast_1d(obs_save)
            obs_save = obs_save.ravel()
            print(obs_log.shape, obs_save.shape)
            return np.concatenate((obs_log, obs_save))
        else:
            return obs_log


class CanonicalMonteCarlo:
    def __init__(self, model, kT, config):
        """

        Parameters
        ----------
        model: dft_latgas object
            DFT lattice gas mapping model
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
        if dE < 0.0:
            self.config = self.model.newconfig(self.config, dconfig)
            self.energy += dE
        else:
            accept_probability = exp(-dE / self.kT)
            if random() <= accept_probability:
                self.config = self.model.newconfig(self.config, dconfig)
                self.energy += dE

    def run(
        self,
        nsteps,
        sample_frequency=verylargeint,
        print_frequency=verylargeint,
        observer=observer_base(),
        save_obs=False,
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
        obs_decode: list
        """
        observables = 0.0
        nsample = 0
        self.energy = self.model.energy(self.config)
        output = open("obs.dat", "a")
        with open(os.devnull, 'w') as f:
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
            args_info = observer.obs_info(self)
            return obs_decode(args_info, observables)
        else:
            return None


# def swap_configs(MCreplicas, rep, accept_count):
#     """
#
#     Parameters
#     ----------
#     MCreplicas:
#     rep:
#     accept_count:
#
#     Returns
#     -------
#
#     """
#     # swap configs, energy
#     tmp = MCreplicas[rep + 1].config
#     tmpe = MCreplicas[rep + 1].energy
#     MCreplicas[rep + 1].config = MCreplicas[rep].config
#     MCreplicas[rep + 1].energy = MCreplicas[rep].energy
#     MCreplicas[rep].config = tmp
#     MCreplicas[rep].energy = tmpe
#     accept_count += 1
#     return MCreplicas, accept_count


"""
class TemperatureReplicaExchange:

    def __init__(self, model, kTs, configs, MCalgo, swap_algo=swap_configs):
        assert len(kTs) == len(configs)
        self.model = model
        self.kTs = kTs
        self.betas = 1.0/kTs
        self.n_replicas = len(kTs)
        self.MCreplicas = []
        self.accept_count = 0
        self.swap_algo = swap_algo
        self.writefunc = writefunc
        self.configs = configs
        for i in range(self.n_replicas):
            self.MCreplicas.append(MCalgo(model, kTs[i], configs[i], writefunc))

    def Xtrial(self):
        # pick a replica
        rep = randrange(self.n_replicas - 1)
        delta = (self.betas[rep + 1] - self.betas[rep]) \
                *(self.MCreplicas[rep].energy -
                  self.MCreplicas[rep+1].energy)
        #print self.MCreplicas[rep].energy, self.model.energy(self.MCreplicas[rep].config)

        if delta < 0.0:
            self.MCreplicas, self.accept_count = self.swap_algo(self.MCreplicas,
                                                           rep, self.accept_count)
            print("RXtrial accepted")
        else:
            accept_probability = exp(-delta)
            #print accept_probability, "accept prob"
            if random() <= accept_probability:
                self.MCreplicas, self.accept_count = self.swap_algo(self.MCreplicas,
                                                           rep, self.accept_count)
                print("RXtrial accepted")
            else:
                print("RXtrial rejected")

    def run(self, nsteps, RXtrial_frequency, pool, sample_frequency=0, subdirs=False):
        self.accept_count = 0
        outerloop = nsteps//RXtrial_frequency
        for i in range(outerloop):
            self.MCreplicas = MultiProcessReplicaRun(self.MCreplicas, RXtrial_frequency, pool, sample_frequency, subdirs)
            self.Xtrial()
        self.configs = [MCreplica.config for MCreplica in self.MCreplicas]
        #print self.accept_count
        #self.accept_count = 0
"""
