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
        config

        Returns
        -------

        """
        return 0.0

    def trialstep(self, config, energy):
        """
        Define a trial step on config. Returns dconfig, which can contain the minimal information for
        constructing the trial configuration from config to be used in newconfig(). Make sure that
        config is the same upon entry and exit

        Parameters
        ----------
        config
        energy

        Returns
        -------

        """
        dE = 0.0
        dconfig = None
        # Return only change in configuration dconfig so that
        # you don't have to copy entire configurations,
        # which can sometimes be costly
        return dconfig, dE

    def newconfig(self, config, dconfig):
        """
        Build new configuration from config and dconfig

        Parameters
        ----------
        config
        dconfig

        Returns
        -------

        """
        return config


class grid_1D:
    def __init__(self, dx, minx, maxx):
        """

        Parameters
        ----------
        dx
        minx
        maxx
        """
        self.dx = dx
        self.x = np.arange(minx, maxx, dx)


def binning(x, nlevels):
    """

    Parameters
    ----------
    x
    nlevels

    Returns
    -------

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
    args

    Returns
    -------

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
    args

    Returns
    -------

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
    args_info
    obs_array

    Returns
    -------

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


def make_observefunc(logfunc, *multiDfuncs):
    def observefunc(calc_state, outputfi):
        """

        Parameters
        ----------
        calc_state
        outputfi

        Returns
        -------

        """
        obs_log = logfunc(calc_state)
        outputfi.write(str(calc_state.kT) + "\t")
        if hasattr(obs_log, "__getitem__"):
            outputfi.write(
                "\t".join([str(observable) for observable in obs_log]) + "\n"
            )
        else:
            outputfi.write(str(obs_log) + "\n")
        obs_ND = []
        for func in multiDfuncs:
            obs_ND.append(func(calc_state))
        return obs_encode(*obs_log, *obs_ND)

    return observefunc


class observer_base:
    def __init__(self):
        self.lprintcount = 0

    def obs_info(self, calc_state):
        """

        Parameters
        ----------
        calc_state

        Returns
        -------

        """
        obs_log = self.logfunc(calc_state)
        if isinstance(obs_log, tuple) == False:
            obs_log = (obs_log,)
        obs_ND = []
        obs_save = self.savefuncs(calc_state)
        if obs_save != None:
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
        calc_state

        Returns
        -------

        """
        return (calc_state.energy,)

    def savefuncs(self, calc_state):
        """

        Parameters
        ----------
        calc_state

        Returns
        -------

        """
        return None

    def writefile(self, calc_state):
        """

        Parameters
        ----------
        calc_state

        Returns
        -------

        """
        return None

    def observe(self, calc_state, outputfi, lprint=True):
        """

        Parameters
        ----------
        calc_state
        outputfi
        lprint

        Returns
        -------

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
        if obs_save != None:
            obs_save = np.atleast_1d(obs_save)
            obs_save = obs_save.ravel()
            print(obs_log.shape, obs_save.shape)
            return np.concatenate((obs_log, obs_save))
        else:
            return obs_log


class CanonicalMonteCarlo:
    def __init__(self, model, kT, config, grid=0):
        """

        Parameters
        ----------
        model
        kT
        config
        grid
        """
        self.model = model
        self.config = config
        self.kT = kT
        self.grid = grid
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
            # print "trial accepted"
        else:
            accept_probability = exp(-dE / self.kT)
            # if sfmt_random.random_sample() <= accept_probability:
            if random() <= accept_probability:
                self.config = self.model.newconfig(self.config, dconfig)
                self.energy += dE
                # print "trial accepted"

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
        nsteps
        sample_frequency
        print_frequency
        observer
        save_obs

        Returns
        -------

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
            if i % sample_frequency == 0 and observe:
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


def swap_configs(MCreplicas, rep, accept_count):
    """

    Parameters
    ----------
    MCreplicas
    rep
    accept_count

    Returns
    -------

    """
    # swap configs, energy
    tmp = MCreplicas[rep + 1].config
    tmpe = MCreplicas[rep + 1].energy
    MCreplicas[rep + 1].config = MCreplicas[rep].config
    MCreplicas[rep + 1].energy = MCreplicas[rep].energy
    MCreplicas[rep].config = tmp
    MCreplicas[rep].energy = tmpe
    accept_count += 1
    return MCreplicas, accept_count


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
