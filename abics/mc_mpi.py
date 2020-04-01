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

import os
import random as rand
from shutil import move
import sys

from mpi4py import MPI

import numpy as np

from abics.mc import observer_base, obs_decode, verylargeint
from abics.util import pickle_dump, pickle_load, numpy_save, numpy_load


class RXParams:
    """ Parameter set for replica exchange Monte Carlo

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
    RXtrial_frequency :
        The number of MC steps between replica exchange operations
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
        self.kTstart = None
        self.kTend = None
        self.nsteps = None
        self.RXtrial_frequency = 1
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
               Dictionary including parameters for replica exchange Monte Carlo method

           Returns
           -------
           params: DFTParams object
               self
        """
        if 'replica' in d:
            d = d['replica']
        params = cls()
        params.nreplicas = d["nreplicas"]
        params.nprocs_per_replica = d["nprocs_per_replica"]
        params.kTstart = d["kTstart"]
        params.kTend = d["kTend"]
        params.nsteps = d["nsteps"]
        params.RXtrial_frequency = d.get("RXtrial_frequency", 1)
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


def RX_MPI_init(rxparams):
    """

    Parameters
    ----------
    rxparams: RXParams
        Parameters for replica exchange Monte Carlo method.

    Returns:
    -------
    comm: comm world
        MPI communicator
    """

    nreplicas = rxparams.nreplicas
    commworld = MPI.COMM_WORLD
    worldrank = commworld.Get_rank()
    worldprocs = commworld.Get_size()
    if rxparams.seed > 0:
        rand.seed(rxparams.seed + worldrank * 137)
    else:
        rand_seeds = [rand.random() for i in range(worldprocs)]
        rand.seed(rand_seeds[worldrank])

    if worldprocs > nreplicas:
        if worldrank == 0:
            print(
                "Setting number of replicas smaller than MPI processes; I hope you"
                + " know what you're doing..."
            )
            sys.stdout.flush()
        if worldrank >= nreplicas:
            # belong to comm that does nothing
            comm = commworld.Split(color=1, key=worldrank)
            comm.Free()
            sys.exit()  # Wait for MPI_finalize
        else:
            comm = commworld.Split(color=0, key=worldrank)
    else:
        comm = commworld
    return comm


class ParallelMC(object):
    def __init__(self, comm, MCalgo, model, configs, kTs, subdirs=True):
        """

        Parameters
        ----------
        comm: comm world
            MPI communicator
        MCalgo: object for MonteCarlo algorithm
            MonteCarlo algorithm
        model: dft_latgas object
            DFT lattice gas mapping model
        configs: config object
            Configurations
        kTs: list
            Temperature list
        subdirs: boolean
            if true,  working directory for this rank is made
        """
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.procs = self.comm.Get_size()
        self.kTs = kTs
        self.model = model
        self.subdirs = subdirs
        self.nreplicas = len(configs)

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

    def run(self, nsteps, sample_frequency, observer=observer_base()):
        """

        Parameters
        ----------
        nsteps: int
            Number of Monte Carlo steps for running.
        sample_frequency: int
            Number of Monte Carlo steps for running.
        observer: observer object

        Returns
        -------
        obs_buffer: numpy array
            Observables
        """
        if self.subdirs:
            # make working directory for this rank
            try:
                os.mkdir(str(self.rank))
            except FileExistsError:
                pass
            os.chdir(str(self.rank))
        observables = self.mycalc.run(nsteps, sample_frequency, observer)
        pickle_dump(self.mycalc.config, "config.pickle")
        if self.subdirs:
            os.chdir("../")
        if sample_frequency:
            obs_buffer = np.empty([self.procs, len(observables)])
            self.comm.Allgather(observables, obs_buffer)
            return obs_buffer


class TemperatureRX_MPI(ParallelMC):
    def __init__(self, comm, MCalgo, model, configs, kTs, subdirs=True):
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
        kTs: list
            Temperature list
        subdirs: boolean
            If true, working directory for this rank is made
        """
        super(TemperatureRX_MPI, self).__init__(
            comm, MCalgo, model, configs, kTs, subdirs
        )
        self.betas = 1.0 / np.array(kTs)
        self.rank_to_T = np.arange(0, self.procs, 1, dtype=np.int)
        self.float_buffer = np.array(0.0, dtype=np.float)
        self.int_buffer = np.array(0, dtype=np.int)
        self.obs_save = []
        self.Trank_hist = []
        self.kT_hist = []

    def reload(self):
        self.rank_to_T = pickle_load("rank_to_T.pickle")
        self.mycalc.kT = self.kTs[self.rank_to_T[self.rank]]
        self.mycalc.config = pickle_load(os.path.join(str(self.rank), "calc.pickle"))
        self.obs_save0 = numpy_load(os.path.join(str(self.rank), "obs_save.npy"))
        self.Trank_hist0 = numpy_load(os.path.join(str(self.rank), "Trank_hist.npy"))
        self.kT_hist0 = numpy_load(os.path.join(str(self.rank), "kT_hist.npy"))
        rand_state = pickle_load(os.path.join(str(self.rank), "rand_state.pickle"))
        rand.setstate(rand_state)

    def find_procrank_from_Trank(self, Trank):
        """

        Parameters
        ----------
        Trank: int
            Temperature rank
        Returns
        -------
        procrank: int
        """
        i = np.argwhere(self.rank_to_T == Trank)
        if i is None:
            sys.exit("Internal error in TemperatureRX_MPI.find_procrank_from_Trank")
        else:
            return i

    def Xtrial(self, XCscheme=-1):
        """

        Parameters
        ----------
        XCscheme: int
        Returns:
        -------

        """
        # What is my temperature rank?
        myTrank = self.rank_to_T[self.rank]
        if (myTrank + XCscheme) % 2 == 0 and myTrank == self.procs - 1:
            self.comm.Allgather(self.rank_to_T[self.rank], self.rank_to_T)
            return
        if XCscheme == 1 and myTrank == 0:
            self.comm.Allgather(self.rank_to_T[self.rank], self.rank_to_T)
            return
        if (myTrank + XCscheme) % 2 == 0:
            myTrankp1 = myTrank + 1
            # Get the energy from the replica with higher temperature
            exchange_rank = self.find_procrank_from_Trank(myTrankp1)
            self.comm.Recv(self.float_buffer, source=exchange_rank, tag=1)
            delta = (self.betas[myTrankp1] - self.betas[myTrank]) * (
                self.mycalc.energy - self.float_buffer
            )
            if delta < 0.0:
                # Exchange temperatures!
                self.comm.Send(
                    [self.rank_to_T[self.rank], 1, MPI.INT], dest=exchange_rank, tag=2
                )

                self.rank_to_T[self.rank] = myTrankp1
            else:
                accept_probability = np.exp(-delta)
                # print accept_probability, "accept prob"
                if rand.random() <= accept_probability:
                    self.comm.Send(
                        [self.rank_to_T[self.rank], 1, MPI.INT],
                        dest=exchange_rank,
                        tag=2,
                    )
                    self.rank_to_T[self.rank] = myTrankp1
                else:
                    # print("RXtrial rejected")
                    self.comm.Send(
                        [self.rank_to_T[exchange_rank], 1, MPI.INT],
                        dest=exchange_rank,
                        tag=2,
                    )
        else:
            myTrankm1 = myTrank - 1
            exchange_rank = self.find_procrank_from_Trank(myTrankm1)
            self.comm.Send(self.mycalc.energy, dest=exchange_rank, tag=1)
            self.comm.Recv([self.int_buffer, 1, MPI.INT], source=exchange_rank, tag=2)
            self.rank_to_T[self.rank] = self.int_buffer
        self.comm.Allgather(self.rank_to_T[self.rank], self.rank_to_T)
        self.mycalc.kT = self.kTs[self.rank_to_T[self.rank]]
        return

    def run(
        self,
        nsteps,
        RXtrial_frequency,
        sample_frequency=verylargeint,
        print_frequency=verylargeint,
        observer=observer_base(),
        subdirs=True,
        save_obs=True,
    ):
        """

        Parameters
        ----------
        nsteps: int
            The number of Monte Carlo steps for running.
        RXtrial_frequency: int
            The number of Monte Carlo steps for replica exchange.
        sample_frequency: int
            The number of Monte Carlo steps for observation of physical quantities.
        print_frequency: int
            The number of Monte Carlo steps for saving physical quantities.
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
        XCscheme = 0
        with open("obs.dat", "a") as output:
            for i in range(1, nsteps + 1):
                self.mycalc.MCstep()
                if i % RXtrial_frequency == 0:
                    self.Xtrial(XCscheme)
                    XCscheme = (XCscheme + 1) % 2
                if observe and i % sample_frequency == 0:
                    obs_step = observer.observe(
                        self.mycalc, output, i % print_frequency == 0
                    )
                    obs[self.rank_to_T[self.rank]] += obs_step
                    if save_obs:
                        self.obs_save.append(obs_step)
                        self.Trank_hist.append(self.rank_to_T[self.rank])
                        self.kT_hist.append(self.mycalc.kT)
                    nsample += 1

                self.comm.Barrier()

                # save information for restart
                pickle_dump(self.mycalc.config, "calc.pickle")
                rand_state = rand.getstate()
                pickle_dump(rand_state, "rand_state.pickle")
                if save_obs:
                    if hasattr(self, "obs_save0"):
                        obs_save_ = np.concatenate((self.obs_save0, np.array(self.obs_save)))
                        Trank_hist_ = np.concatenate(
                            (self.Trank_hist0, np.array(self.Trank_hist))
                        )
                        kT_hist_ = np.concatenate((self.kT_hist0, np.array(self.kT_hist)))
                    else:
                        obs_save_ = np.array(self.obs_save)
                        Trank_hist_ = np.array(self.Trank_hist)
                        kT_hist_ = np.array(self.kT_hist)

                    numpy_save(obs_save_, "obs_save.npy")
                    numpy_save(Trank_hist_, "Trank_hist.npy")
                    numpy_save(kT_hist_, "kT_hist.npy")

                if subdirs:
                    os.chdir("../")
                if self.rank == 0:
                    pickle_dump(self.rank_to_T, "rank_to_T.pickle")
                    numpy_save(self.kTs, "kTs.npy")
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
            args_info = observer.obs_info(self.mycalc)
            for i in range(len(self.kTs)):
                obs_list.append(obs_decode(args_info, obs_buffer[i]))
            return obs_list
