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

from typing import Type

import os
import sys

from mpi4py import MPI

import numpy as np
import numpy.random as rand

from abics.mc import ObserverBase, verylargeint, MCAlgorithm, Model
from abics.util import pickle_dump, pickle_load, numpy_save, numpy_load


class SamplerParams:
    """Parameter set for specifying sampling algorithm

    Attributes
    ----------
    sampler : str
        Sampler name
    """

    def __init__(self):
        self.sampler = "RXMC"

    @classmethod
    def from_dict(cls, d):
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary including parameters for specifying sampling algorithm

        Returns
        -------
        params: SamplerParams object
            self
        """
        if "sampler" in d:
            d = d["sampler"]
        params = cls()
        params.sampler = d.get("sampler", "RXMC")
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
        SamplerParams: SamplerParams object
            self
        """
        import toml

        return cls.from_dict(toml.load(fname))


class RefParams:
    """Parameter set for reference calculations in active learning run

    Attributes
    ----------
    nreplicas : int
        The number of replicas
    ndata: int
        The number of structures to be converted
    """

    nreplicas: int
    "The number of replicas"
    ndata: int
    "The number of structures to be converted"
    sampler: str
    "The sampling method ('linspace' or 'random')"

    def __init__(self):
        self.nreplicas = 1
        self.ndata = 1
        self.sampler = "linspace"
        self.seed = 0

    def sampling(self, nsamples: int) -> np.ndarray:
        if self.sampler == "linspace":
            ret = np.linspace(0, nsamples - 1, num=self.ndata, dtype=int)
            return ret
        elif self.sampler == "random":
            ret = rand.choice(np.arange(nsamples), size=self.ndata, replace=False)
            ret.sort()
            return ret
        else:
            print(f"Unknown sampler: {self.sampler}")
            sys.exit(1)

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
        params = cls()
        params.nreplicas = d["nreplicas"]
        params.ndata = d["ndata"]
        params.sampler = d.get("sampler", "linspace")
        # params.nprocs_per_replica = d["nprocs_per_replica"]
        # params.nsteps = d["nsteps"]
        # params.sample_frequency = d.get("sample_frequency", 1)
        # params.reload = d.get("reload", False)
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

        d = toml.load(fname)
        return cls.from_dict(d["mlref"])


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
        self.nsteps = None
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


class ParalleMCParams:
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
        self.kTstart = None
        self.kTend = None
        self.nsteps = None
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


class RXParams:
    """Parameter set for replica exchange Monte Carlo

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

        d = toml.load(fname)
        return cls.from_dict(d["sampling"])


def RX_MPI_init(rxparams: RXParams, dftparams=None):
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
    if dftparams != None and dftparams.ensemble and dftparams.par_ensemble:
        nensemble = len(dftparams.base_input_dir)
    else:
        nensemble = 1
    nreplicas = rxparams.nreplicas * nensemble
    commworld = MPI.COMM_WORLD
    worldrank = commworld.Get_rank()
    worldprocs = commworld.Get_size()

    if worldprocs < nreplicas:
        if worldrank == 0:
            print(
                "ERROR! Please run with at least as many MPI processes as the number of replicas"
            )
        sys.exit(1)

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
    comm = comm.Create_cart(
        dims=[rxparams.nreplicas, nensemble], periods=[False, False], reorder=True
    )
    commRX = comm.Sub(remain_dims=[True, False])
    commEnsemble = comm.Sub(remain_dims=[False, True])
    RXrank = commRX.Get_rank()
    if rxparams.seed > 0:
        rand.seed(rxparams.seed + RXrank * 137)
    else:
        rand_seeds = [rand.randint(10000) for i in range(commRX.Get_size())]
        rand_seed = commEnsemble.bcast(rand_seeds[RXrank], root=0)
        rand.seed(rand_seed)

    # return commRX
    if dftparams == None:
        return commRX
    return commRX, commEnsemble, comm


class ParallelMC(object):
    def __init__(
        self,
        comm,
        MCalgo: Type[MCAlgorithm],
        model: Model,
        configs,
        kTs,
        subdirs=True,
        write_node=True,
    ):
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

    def run(self, nsteps, sample_frequency, observer=ObserverBase()):
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
        if self.write_node:
            pickle_dump(self.mycalc.config, "config.pickle")
        if self.subdirs:
            os.chdir("../")
        if sample_frequency:
            obs_buffer = np.empty([self.procs, len(observables)])
            self.comm.Allgather(observables, obs_buffer)
            return obs_buffer


class EmbarrassinglyParallelSampling:
    def __init__(
        self,
        comm,
        MCalgo: Type[MCAlgorithm],
        model: Model,
        configs,
        kTs=None,
        subdirs=True,
        write_node=True,
    ):
        """

        Parameters
        ----------
        comm: comm world
            MPI communicator
        MCalgo: Type[MCAlgorithm]
            MonteCarlo algorithm class (not instance)
        model: Model
            Model
        configs: config object
            Configuration
        kTs: list
            Temperature list
        subdirs: boolean
            If true, working directory for this rank is made
        """
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.procs = self.comm.Get_size()
        if kTs is None:
            kTs = [0] * self.procs
        self.kTs = kTs
        self.model = model
        self.subdirs = subdirs
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
        self.obs_save = []
        self.kT_hist = []
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
        nsteps,
        sample_frequency=verylargeint,
        print_frequency=verylargeint,
        observer=ObserverBase(),
        subdirs=True,
        save_obs=True,
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
            for i in range(1, nsteps + 1):
                self.mycalc.MCstep()

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
    def __init__(self, comm, MCalgo, model, configs, subdirs=True, write_node=True):
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
        self.write_node = write_node


class TemperatureRX_MPI(ParallelMC):
    def __init__(
        self, comm, MCalgo, model, configs, kTs, subdirs=True, write_node=True
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
        kTs: list
            Temperature list
        subdirs: boolean
            If true, working directory for this rank is made
        """
        super(TemperatureRX_MPI, self).__init__(
            comm, MCalgo, model, configs, kTs, subdirs, write_node
        )
        self.betas = 1.0 / np.array(kTs)
        self.rank_to_T = np.arange(0, self.procs, 1, dtype=np.int)
        self.float_buffer = np.array(0.0, dtype=np.float)
        self.int_buffer = np.array(0, dtype=np.int)
        self.obs_save = []
        self.Trank_hist = []
        self.kT_hist = []
        self.Lreload = False

    def reload(self):
        self.rank_to_T = pickle_load("rank_to_T.pickle")
        self.mycalc.kT = self.kTs[self.rank_to_T[self.rank]]
        self.mycalc.config = pickle_load(os.path.join(str(self.rank), "calc.pickle"))
        self.obs_save0 = numpy_load(os.path.join(str(self.rank), "obs_save.npy"))
        self.mycalc.energy = self.obs_save0[-1, 0]
        self.Trank_hist0 = numpy_load(os.path.join(str(self.rank), "Trank_hist.npy"))
        self.kT_hist0 = numpy_load(os.path.join(str(self.rank), "kT_hist.npy"))
        rand_state = pickle_load(os.path.join(str(self.rank), "rand_state.pickle"))
        rand.set_state(rand_state)
        self.Lreload = True

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
        observer=ObserverBase(),
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
        XCscheme = 0
        with open("obs.dat", "a") as output:
            for i in range(1, nsteps + 1):
                self.mycalc.MCstep()
                if i % RXtrial_frequency == 0:
                    self.Xtrial(XCscheme)
                    XCscheme = (XCscheme + 1) % 2
                if observe and i % sample_frequency == 0:
                    obs_step = observer.observe(
                        self.mycalc,
                        output,
                        i % print_frequency == 0 and self.write_node,
                    )
                    obs[self.rank_to_T[self.rank]] += obs_step
                    if save_obs:
                        self.obs_save.append(obs_step)
                        self.Trank_hist.append(self.rank_to_T[self.rank])
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
                                Trank_hist_ = np.concatenate(
                                    (self.Trank_hist0, np.array(self.Trank_hist))
                                )
                                kT_hist_ = np.concatenate(
                                    (self.kT_hist0, np.array(self.kT_hist))
                                )
                            else:
                                obs_save_ = np.array(self.obs_save)
                                Trank_hist_ = np.array(self.Trank_hist)
                                kT_hist_ = np.array(self.kT_hist)

                            numpy_save(obs_save_, "obs_save.npy")
                            numpy_save(Trank_hist_, "Trank_hist.npy")
                            numpy_save(kT_hist_, "kT_hist.npy")

                    if subdirs:
                        os.chdir("../")
                    if self.write_node:
                        if self.rank == 0:
                            pickle_dump(self.rank_to_T, "rank_to_T.pickle")
                            numpy_save(self.kTs, "kTs.npy")
                    if subdirs:
                        os.chdir(str(self.rank))

        if nsample != 0:
            obs = np.array(obs)
            obs_buffer = np.empty(obs.shape)
            obs /= nsample
            self.comm.Allreduce(obs, obs_buffer, op=MPI.SUM)
            obs_list = []
            obs_info = observer.obs_info(self.mycalc)
            for i in range(len(self.kTs)):
                obs_list.append(obs_info.decode(obs_buffer[i]))
            if subdirs:
                os.chdir("../")
            return obs_list

        if subdirs:
            os.chdir("../")
