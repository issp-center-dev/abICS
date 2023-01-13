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
from typing import overload

import os
import sys

from mpi4py import MPI

import numpy as np
import numpy.random as rand

from abics.observer import ObserverBase
from abics.model import Model
from abics.sampling.mc import verylargeint, MCAlgorithm
from abics.util import pickle_dump


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
            ret = np.linspace(0, nsamples - 1, num=self.ndata, dtype=np.int64)
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


@overload
def RX_MPI_init(nreplicas: int) -> MPI.Cartcomm:
    ...


@overload
def RX_MPI_init(nreplicas: int, seed: int) -> MPI.Cartcomm:
    ...


@overload
def RX_MPI_init(nreplicas: int, seed: int, nensemble: None) -> MPI.Cartcomm:
    ...


@overload
def RX_MPI_init(
    nreplicas: int, seed: int, nensemble: int
) -> tuple[MPI.Cartcomm, MPI.Cartcomm, MPI.Cartcomm]:
    ...


def RX_MPI_init(nreplicas: int, seed: int = 0, nensemble=None):
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
    if nensemble is None:
        nensemble_ = 1
    else:
        nensemble_ = nensemble
    nreplicas = nreplicas * nensemble_
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
            comm = MPI.Intracomm(commworld.Split(color=0, key=worldrank))
            # comm = commworld.Split(color=0, key=worldrank)
    else:
        comm = commworld
    comm = comm.Create_cart(
        dims=[nreplicas, nensemble_], periods=[False, False], reorder=True
    )
    commRX = comm.Sub(remain_dims=[True, False])
    commEnsemble = comm.Sub(remain_dims=[False, True])
    RXrank = commRX.Get_rank()
    if seed > 0:
        rand.seed(seed + RXrank * 137)
    else:
        rand_seeds = [rand.randint(10000) for i in range(commRX.Get_size())]
        rand_seed = commEnsemble.bcast(rand_seeds[RXrank], root=0)
        rand.seed(rand_seed)

    # return commRX
    if nensemble is None:
        return commRX
    else:
        return commRX, commEnsemble, comm


class ParallelMC(object):
    comm: MPI.Comm
    "MPI Communicator"
    rank: int
    "MPI rank"
    procs: int
    "MPI size"
    kTs: list[float]
    "Temperatures"
    nreplicas: int
    "The number of replicas"
    model: Model
    mycalc: MCAlgorithm
    write_node: bool
    obs_save: list[np.ndarray]
    "Observed data (e.g., energy)"
    obs_save0: np.ndarray
    "Observed data before reloading"
    Lreload: bool
    "Have reloaded or not"

    def __init__(
        self,
        comm,
        MCalgo: type[MCAlgorithm],
        model: Model,
        configs,
        kTs,
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
        """
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.procs = self.comm.Get_size()
        self.kTs = kTs
        self.model = model
        self.nreplicas = len(configs)
        self.write_node = write_node

        ## mycalc.kT and mycalc.config should be set later
        myconfig = configs[0]
        mytemp = kTs[0]
        self.mycalc = MCalgo(model, mytemp, myconfig)

    def run(
        self,
        nsteps: int,
        sample_frequency: int = verylargeint,
        nsubsteps_in_step: int = 1,
        observer: ObserverBase = ObserverBase(),
        subdirs: bool = True,
    ):
        """

        Parameters
        ----------
        nsteps: int
            Number of Monte Carlo steps for running.
        sample_frequency: int
            Number of Monte Carlo steps for running.
        nsubsteps_in_step: int
            Number of Monte Carlo substeps in one MC step.
        observer: observer object
        subdirs: boolean
            if true,  working directory for this rank is made

        Returns
        -------
        obs_buffer: numpy array
            Observables
        """
        if subdirs:
            # make working directory for this rank
            try:
                os.mkdir(str(self.rank))
            except FileExistsError:
                pass
            os.chdir(str(self.rank))
        observables = self.mycalc.run(
            nsteps,
            sample_frequency=sample_frequency,
            nsubsteps_in_step=nsubsteps_in_step,
            observer=observer,
        )
        if self.write_node:
            pickle_dump(self.mycalc.config, "config.pickle")
        if subdirs:
            os.chdir("../")
        if sample_frequency:
            obs_buffer = np.empty([self.procs, len(observables)])
            self.comm.Allgather(observables, obs_buffer)
            return obs_buffer
