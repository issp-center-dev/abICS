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

import os
import sys

from mpi4py import MPI

import numpy as np
import numpy.random as rand

from abics.model import Model
from abics.observer import ObserverBase
from abics.sampling.mc import MCAlgorithm, verylargeint, write_obs_header
from abics.sampling.mc_mpi import ParallelMC
from abics.util import pickle_dump, pickle_load, numpy_save, numpy_load


class ParallelRandomParams:
    """Parameter set for parallel random sampling

    Attributes
    ----------
    nreplicas : int
        The number of replicas
    nprocs_per_replica : int
        The number of processes which a replica uses
    nsteps : int
        The number of MC steps
    sample_frequency :
        The number of MC steps between measurements observables
    print_frequency :
        The number of MC steps between show information
    reload : bool
        Whether to restart simulation or not
    seed : int
        The seed of the random number generator
        If 0, some random number is used (e.g., system time or some random noise).
    """

    def __init__(self):
        self.nreplicas = None
        self.nprocs_per_replica = 1
        self.nsteps = 0
        self.sample_frequency = 1
        self.print_frequency = 1
        self.reload = False
        self.seed = 0

    @classmethod
    def from_dict(cls, d):
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary including parameters for parallel random sampling

        Returns
        -------
        params: DFTParams object
            self
        """
        if "replica" in d:
            d = d["replica"]
        params = cls()
        params.nreplicas = d["nreplicas"]
        params.nprocs_per_replica = d["nprocs_per_replica"]
        params.nsteps = d["nsteps"]
        params.sample_frequency = d.get("sample_frequency", 1)
        params.print_frequency = d.get("print_frequency", 1)
        params.reload = d.get("reload", False)
        params.seed = d.get("seed", 0)
        return params

    @classmethod
    def from_toml(cls, fname):
        """
        Read information from toml file

        Parameters
        ----------
        f: str
            The name of input toml File

        Returns
        -------
        DFTParams: DFTParams object
            self
        """
        import toml

        return cls.from_dict(toml.load(fname))


class ParallelMCParams:
    """Parameter set for embarrasingly parallel Monte Carlo

    Attributes
    ----------
    nreplicas : int
        The number of replicas
    nprocs_per_replica : int
        The number of processes which a replica uses
    kTstart : float
        The lower bound of temperature range
    kTend : float
        The upper bound of temperature range
    nsteps : int
        The number of MC steps
    sample_frequency :
        The number of MC steps between measurements observables
    print_frequency :
        The number of MC steps between show information
    reload : bool
        Whether to restart simulation or not
    seed : int
        The seed of the random number generator
        If 0, some random number is used (e.g., system time or some random noise).
    """

    def __init__(self):
        self.nreplicas = None
        self.nprocs_per_replica = 1
        self.kTstart = 0.0
        self.kTend = 1.0
        self.nsteps = 0
        self.sample_frequency = 1
        self.print_frequency = 1
        self.reload = False
        self.seed = 0

    @classmethod
    def from_dict(cls, d):
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary including parameters for embarrassingly parallel Monte Carlo method

        Returns
        -------
        params: DFTParams object
            self
        """
        if "replica" in d:
            d = d["replica"]
        params = cls()
        params.nreplicas = d["nreplicas"]
        params.nprocs_per_replica = d["nprocs_per_replica"]
        params.kTstart = d["kTstart"]
        params.kTend = d["kTend"]
        params.nsteps = d["nsteps"]
        params.sample_frequency = d.get("sample_frequency", 1)
        params.print_frequency = d.get("print_frequency", 1)
        params.reload = d.get("reload", False)
        params.seed = d.get("seed", 0)
        return params

    @classmethod
    def from_toml(cls, fname):
        """
        Read information from toml file

        Parameters
        ----------
        f: str
            The name of input toml File

        Returns
        -------
        DFTParams: DFTParams object
            self
        """
        import toml

        return cls.from_dict(toml.load(fname))


class EmbarrassinglyParallelSampling:
    def __init__(
        self,
        comm,
        MCalgo: type[MCAlgorithm],
        model: Model,
        configs,
        kTs=None,
        write_node=True,
    ):
        """

        Parameters
        ----------
        comm: comm world
            MPI communicator
        MCalgo: type[MCAlgorithm]
            MonteCarlo algorithm class (not instance)
        model: Model
            Model
        configs: config object
            Configuration
        kTs: list
            Temperature list
        """
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.procs = self.comm.Get_size()
        if kTs is None:
            kTs = [0] * self.procs
        self.kTs = kTs
        self.model = model
        self.nreplicas = len(configs)
        self.write_node = write_node

        if not (self.procs == self.nreplicas == len(self.kTs)):
            if self.rank == 0:
                print(
                    "ERROR: You have to set the number of replicas equal to the"
                    + "number of temperatures equal to the number of processes"
                )
            sys.exit(1)

        myconfig = configs[self.rank]
        mytemp = kTs[self.rank]
        self.mycalc = MCalgo(model, mytemp, myconfig)
        self.obs_save: list[np.ndarray] = []
        self.kT_hist: list[float] = []
        self.Lreload = False

    def reload(self):
        self.mycalc.config = pickle_load(os.path.join(str(self.rank), "calc.pickle"))
        self.obs_save0 = numpy_load(os.path.join(str(self.rank), "obs_save.npy"))
        self.mycalc.energy = self.obs_save0[-1, 0]
        self.kT_hist0 = numpy_load(os.path.join(str(self.rank), "kT_hist.npy"))
        self.mycalc.kT = self.kT_hist0[-1]
        rand_state = pickle_load(os.path.join(str(self.rank), "rand_state.pickle"))
        rand.set_state(rand_state)
        self.Lreload = True

    def run(
        self,
        nsteps: int,
        sample_frequency: int = verylargeint,
        print_frequency: int = verylargeint,
        nsubsteps_in_step: int = 1,
        observer: ObserverBase = ObserverBase(),
        subdirs: bool = True,
        save_obs: bool = True,
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
        subdirs: boolean
            If true, working directory for this rank is made
        save_obs: boolean

        Returns
        -------
        obs_list: list
            Observation list
        """
        if subdirs:
            try:
                os.mkdir(str(self.rank))
            except FileExistsError:
                pass
            os.chdir(str(self.rank))
        self.accept_count = 0
        if not self.Lreload:
            self.mycalc.energy = self.mycalc.model.energy(self.mycalc.config)
        with open(os.devnull, "w") as f:
            test_observe = observer.observe(self.mycalc, f, lprint=False)
        if hasattr(test_observe, "__iter__"):
            obs_len = len(test_observe)
            obs = np.zeros([len(self.kTs), obs_len])
        if hasattr(test_observe, "__add__"):
            observe = True
        else:
            observe = False
        nsample = 0
        with open("obs.dat", "a") as output:
            write_obs_header(output, self.mycalc, observer)
            for i in range(1, nsteps + 1):
                self.mycalc.MCstep(nsubsteps_in_step)

                if observe and i % sample_frequency == 0:
                    obs_step = observer.observe(
                        self.mycalc,
                        output,
                        i % print_frequency == 0 and self.write_node,
                    )
                    obs[self.rank] += obs_step
                    if save_obs:
                        self.obs_save.append(obs_step)
                        self.kT_hist.append(self.mycalc.kT)
                    nsample += 1

                self.comm.Barrier()

                if self.write_node:
                    # save information for restart
                    pickle_dump(self.mycalc.config, "calc.pickle")
                    rand_state = rand.get_state()
                    pickle_dump(rand_state, "rand_state.pickle")
                    if save_obs:
                        if hasattr(self, "obs_save0"):
                            obs_save_ = np.concatenate(
                                (self.obs_save0, np.array(self.obs_save))
                            )
                            kT_hist_ = np.concatenate(
                                (self.kT_hist0, np.array(self.kT_hist))
                            )
                        else:
                            obs_save_ = np.array(self.obs_save)
                            kT_hist_ = np.array(self.kT_hist)

                        numpy_save(obs_save_, "obs_save.npy")
                        numpy_save(kT_hist_, "kT_hist.npy")

                if subdirs:
                    os.chdir("../")
                if subdirs:
                    os.chdir(str(self.rank))

        if subdirs:
            os.chdir("../")

        if nsample != 0:
            obs = np.array(obs)
            obs_buffer = np.empty(obs.shape)
            obs /= nsample
            self.comm.Allreduce(obs, obs_buffer, op=MPI.SUM)
            obs_list = []
            obs_info = observer.obs_info(self.mycalc)
            for i in range(len(self.kTs)):
                obs_list.append(obs_info.decode(obs_buffer[i]))
            return obs_list


class RandomSampling_MPI(ParallelMC):
    def __init__(
        self, comm, MCalgo: type[MCAlgorithm], model: Model, configs, write_node=True
    ):
        """

        Parameters
        ----------
        comm: comm world
            MPI communicator
        MCalgo: object for MonteCarlo algorithm
            MonteCarlo algorithm
        model: dft_latgas
            DFT lattice gas mapping  model
        configs: config object
            Configuration
        """

        super().__init__(comm, MCalgo, model, configs, kTs, write_node=write_node)
        self.mycalc.kT = kTs[self.rank]
        self.mycalc.config = configs[self.rank]
        self.betas = 1.0 / np.array(kTs)
        self.rank_to_T = np.arange(0, self.procs, 1, dtype=np.int64)
        self.float_buffer = np.array(0.0, dtype=np.float64)
        self.int_buffer = np.array(0, dtype=np.int64)
        self.obs_save = []
        self.Trank_hist = []
        self.kT_hist = []
        self.write_node = write_node
