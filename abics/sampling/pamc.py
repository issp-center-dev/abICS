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
from abics.sampling.mc import verylargeint, MCAlgorithm
from abics.sampling.mc_mpi import ParallelMC
from abics.sampling.resampling import WalkerTable
from abics.util import pickle_dump, pickle_load, numpy_save, numpy_load


class PAMCParams:
    """Parameter set for population annealing Monte Carlo

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
    kTnum : int
        The number of temperature points
    nsteps : int
        The number of MC steps between annaling
    resample_frequency :
        The number of annealing between resampling
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
        self.nreplicas = 1
        self.nprocs_per_replica = 1
        self.kTstart = 0.0
        self.kTend = 1.0
        self.kTnum = 1
        self.nsteps = 0
        self.resample_frequency = 1
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
        params.kTnum = d["kTnum"]
        if "nsteps_between_anneal" in d:
            params.nsteps = d["nsteps_between_anneal"] * params.kTnum
            if "nsteps" in d:
                msg = f"Error: Both nsteps and nsteps_between_anneal are specified"
                raise RuntimeError(msg)
        else:
            if "nsteps" in d:
                params.nsteps = d["nsteps"]
            else:
                msg = f"Error: Neither nsteps nor nsteps_between_anneal are specified"
                raise RuntimeError(msg)
        params.resample_frequency = d.get("resample_frequency", 1)
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


class PopulationAnnealing(ParallelMC):
    Tindex: int
    logweight: float
    "Logarithm of Neal-Jarzynski weight"
    logweight_history: list[float]
    dlogz: list[float]
    kT_history: list[float]
    acceptance_ratios: list[float]

    def __init__(
        self,
        comm,
        MCalgo: type[MCAlgorithm],
        model: Model,
        configs,
        kTs: np.ndarray,
        write_node=True,
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
        if kTs[0] < kTs[-1]:
            if comm.Get_rank() == 0:
                print("Warning: kTs[0] < kTs[-1]. abICS reverses kTs.")
                sys.stdout.flush()
            kTs = kTs[::-1]
        if kTs[0] != np.inf:
            kTs = np.concatenate((np.array([np.inf]), kTs))
        super().__init__(comm, MCalgo, model, configs, kTs, write_node=write_node)
        self.mycalc.kT = kTs[0]
        self.mycalc.config = configs[self.rank]
        self.betas = 1.0 / np.array(kTs)
        self.float_buffer = np.array(0.0, dtype=np.float64)
        self.int_buffer = np.array(0, dtype=np.int64)
        self.obs_save = []
        self.Tindex = 0
        self.logweight = 0.0
        self.logweight_history = []
        self.dlogz = []
        self.kT_history = []
        self.acceptance_ratios = []
        self.Lreload = False
        self.use_resample_old = False

    def reload(self):
        self.mycalc.config = pickle_load(os.path.join(str(self.rank), "calc.pickle"))
        self.obs_save0 = numpy_load(os.path.join(str(self.rank), "obs_save.npy"))
        self.mycalc.energy = self.obs_save0[-1, 0]
        wh = numpy_load(os.path.join(str(self.rank), "logweight_hist.npy"))
        self.logweight_history = [w for w in wh]
        self.logweight = self.logweight_history[-1]
        kh = numpy_load(os.path.join(str(self.rank), "kT_hist.npy"))
        self.kT_history = [T for T in kh]
        self.Tindex = np.unique(kh).size
        rand_state = pickle_load(os.path.join(str(self.rank), "rand_state.pickle"))
        rand.set_state(rand_state)
        self.Lreload = True

    def save(self, save_obs: bool):
        # save information for restart
        pickle_dump(self.mycalc.config, "calc.pickle")
        rand_state = rand.get_state()
        pickle_dump(rand_state, "rand_state.pickle")
        if save_obs:
            if hasattr(self, "obs_save0"):
                obs_save_ = np.concatenate((self.obs_save0, np.array(self.obs_save)))
            else:
                obs_save_ = np.array(self.obs_save)
            numpy_save(obs_save_, "obs_save.npy")
            numpy_save(self.logweight_history, "logweight_hist.npy")
            numpy_save(self.kT_history, "kT_hist.npy")

    def anneal(self, energy: float):
        assert 0 < self.Tindex < len(self.kTs)
        mdbeta = self.betas[self.Tindex - 1] - self.betas[self.Tindex]
        dlogz = mdbeta * energy
        self.dlogz.append(dlogz)
        self.logweight += dlogz
        self.mycalc.kT = self.kTs[self.Tindex]

    def resample(self):
        if self.use_resample_old:
            self.__resample_old()
        else:
            self.__resample_new()

    def __resample_old(self):
        buffer = np.zeros((self.procs, 2), dtype=np.float64)
        self.comm.Allgather(np.array((self.logweight, self.mycalc.energy)), buffer)
        logweights = buffer[:, 0]
        energies = buffer[:, 1]

        configs = self.comm.allgather(self.mycalc.config)

        weights = np.exp(logweights - np.max(logweights))  # avoid overflow
        acc_weights = np.add.accumulate(weights)
        r = rand.rand() * acc_weights[-1]
        i = np.searchsorted(acc_weights, r)

        self.mycalc.config = configs[i]
        self.mycalc.energy = energies[i]
        self.logweight = 0.0

    def __resample_new(self):
        buffer = np.zeros((self.procs, 2), dtype=np.float64)
        self.comm.Allgather(np.array((self.logweight, self.mycalc.energy)), buffer)
        logweights = buffer[:, 0]
        energies = buffer[:, 1]

        weights = np.exp(logweights - np.max(logweights))  # avoid overflow
        resampler = WalkerTable(weights)
        index = resampler.sample(size=self.procs)
        self.comm.Bcast(index, root=0)
        mysrc = index[self.rank]
        mydst: list[int] = []
        for dst, src in enumerate(index):
            if src == self.rank:
                mydst.append(dst)
        send_reqs = [
            self.comm.isend(self.mycalc.config, dest=dst, tag=dst) for dst in mydst
        ]
        recv_req = self.comm.irecv(source=mysrc, tag=self.rank)
        for req in send_reqs:
            req.wait()
        newconfig = recv_req.wait()
        self.comm.Barrier()
        self.mycalc.config = newconfig
        self.mycalc.energy = energies[mysrc]
        self.logweight = 0.0

    def run(
        self,
        nsteps: int,
        resample_frequency: int = 1,
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
            The number of Monte Carlo steps.
        resample_frequency: int
            The number of anneals between resampling.
        sample_frequency: int
            The number of Monte Carlo steps for observation of physical quantities.
        print_frequency: int
            The number of Monte Carlo steps for saving physical quantities.
        nsubsteps_in_step: int
            The number of Monte Carlo substeps in one MC step.
        observer: observer object
        subdirs: boolean
            If true, working directory for this rank is made
        save_obs: boolean

        Returns
        -------
        obs_list: list
            Observation list
        """

        kTnum = len(self.kTs)
        nba = nsteps // kTnum
        nsteps_between_anneal = nba * np.ones(kTnum, dtype=int)
        if nsteps - nba*kTnum > 0:
            nsteps_between_anneal[-(nsteps - nba*kTnum):] += 1

        if subdirs:
            try:
                os.mkdir(str(self.rank))
            except FileExistsError:
                pass
            os.chdir(str(self.rank))
        if not self.Lreload:
            self.mycalc.energy = self.mycalc.model.energy(self.mycalc.config)
        with open(os.devnull, "w") as f:
            test_observe = observer.observe(self.mycalc, f, lprint=False)
            obs_len = len(test_observe)
            obs = np.zeros([len(self.kTs), obs_len])
        nsample = 0
        output = open("obs.dat", "a")
        numT = self.betas.size
        while self.Tindex < numT:
            if self.Tindex > 0:
                self.anneal(self.mycalc.energy)
                if self.rank == 0:
                    print(
                        "--Anneal from T={} to {}".format(
                            self.kTs[self.Tindex - 1], self.kTs[self.Tindex]
                        )
                    )
                    sys.stdout.flush()
                if self.Tindex % resample_frequency == 0:
                    self.resample()
                    if self.rank == 0:
                        print("--Resampling finishes")
                        sys.stdout.flush()
            ntrials = self.mycalc.ntrials
            naccepted = self.mycalc.naccepted
            for i in range(1, nsteps_between_anneal[self.Tindex] + 1):
                self.mycalc.MCstep(nsubsteps_in_step)
                if i % sample_frequency == 0:
                    obs_step = observer.observe(
                        self.mycalc,
                        output,
                        i % print_frequency == 0 and self.write_node,
                    )
                    obs[self.Tindex, :] += obs_step
                    if save_obs:
                        self.obs_save.append(obs_step)
                        self.logweight_history.append(self.logweight)
                        self.kT_history.append(self.mycalc.kT)
                    nsample += 1
                    if self.write_node:
                        self.save(save_obs)
            naccepted = self.mycalc.naccepted - naccepted
            ntrials = self.mycalc.ntrials - ntrials
            self.acceptance_ratios.append(naccepted / ntrials)
            with open("acceptance_ratio.dat", "w") as f:
                for T, ar in zip(self.kTs, self.acceptance_ratios):
                    f.write(f"{T} {ar}\n")
            self.Tindex += 1
        output.close()

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
            if save_obs:
                self.postproc()
            return obs_list

        if subdirs:
            os.chdir("../")

    def postproc(self):
        nT = len(self.kTs)
        # logweights = numpy_load(os.path.join(str(self.rank), "logweight_hist.npy"))
        logweights = np.array(self.logweight_history)
        nsamples = logweights.shape[0] // nT
        logweights = np.array(logweights[::nsamples])
        dlogz = np.array(self.dlogz + [0.0])

        if hasattr(self, "obs_save0"):
            obs1 = np.concatenate((self.obs_save0, np.array(self.obs_save)))
        else:
            obs1 = np.array(self.obs_save)
        obs1 = obs1.reshape(nT, nsamples, -1)
        # obs1 = numpy_load(os.path.join(str(self.rank), "obs_save.npy")).reshape(
        # nT, nsamples, -1
        # )
        nobs = obs1.shape[2]
        obs2 = (obs1 * obs1).mean(axis=1)
        obs1 = obs1.mean(axis=1)
        obs = np.zeros((nT, 2 * nobs), dtype=np.float64)
        for iobs in range(nobs):
            obs[:, 2 * iobs] = obs1[:, iobs]
            obs[:, 2 * iobs + 1] = obs2[:, iobs]

        nreplicas = self.procs
        # wall = np.zeros((nreplicas, nT, 2), dtype=np.float64)
        # obs_all = np.zeros((nreplicas, nT, 2 * nobs), dtype=np.float64)
        buffer_all = np.zeros((nreplicas, nT, 2 + 2 * nobs), dtype=np.float64)
        self.comm.Allgather(
            np.concatenate(
                [logweights.reshape(-1, 1), dlogz.reshape(-1, 1), obs], axis=1
            ),
            buffer_all,
        )
        # self.comm.Allgather(obs, obs_all)

        logweights_all = buffer_all[:, :, 0]
        dlogz_all = buffer_all[:, :, 1]
        obs_all = buffer_all[:, :, 2:]

        logweights_max = logweights_all.max(axis=0)
        weights = np.exp(logweights_all - logweights_max)
        ow = np.einsum("ito,it->ito", obs_all, weights)

        lzw = logweights_all + dlogz_all - logweights_max
        lzw_max = lzw.max(axis=0)
        zw = np.exp(lzw - lzw_max)

        # bootstrap
        index = np.random.randint(nreplicas, size=nreplicas)
        numer = ow[index, :, :].mean(axis=0)
        zw_numer = zw[index, :].mean(axis=0)
        denom = weights[index, :].mean(axis=0)
        o = np.zeros((nT, 3 * nobs + 1), dtype=np.float64)
        for iobs in range(nobs):
            o[:, 3 * iobs] = numer[:, 2 * iobs] / denom[:]
            o[:, 3 * iobs + 1] = numer[:, 2 * iobs + 1] / denom[:]
            o[:, 3 * iobs + 2] = o[:, 3 * iobs + 1] - o[:, 3 * iobs] ** 2
        o[:, 3 * nobs] = zw_numer / denom
        o_all = np.zeros((nreplicas, nT, 3 * nobs + 1), dtype=np.float64)
        self.comm.Allgather(o, o_all)
        if self.rank == 0:
            o_mean = o_all.mean(axis=0)
            o_err = o_all.std(axis=0)
            with open("result.dat", "w") as f:
                for iT in range(nT):
                    f.write(str(self.kTs[iT]))
                    for iobs in range(3 * nobs):
                        f.write(f" {o_mean[iT, iobs]} {o_err[iT, iobs]}")
                    f.write("\n")
            dlogZ = np.log(o_mean[:-1, 3 * nobs]) + lzw_max[:-1]
            dlogZ_err = o_err[:-1, 3 * nobs] / o_mean[:-1, 3 * nobs]
            with open("logZ.dat", "w") as f:
                F = 0.0
                dF = 0.0
                f.write(f"{self.kTs[0]} {F} {dF} {0.0} {0.0}\n")
                for i, (dlz, dlz_e) in enumerate(zip(dlogZ, dlogZ_err)):
                    F += dlz
                    dF += dlz_e
                    f.write(f"{self.kTs[i+1]} {F} {dF} {dlz} {dlz_e}\n")

        self.comm.Barrier()
