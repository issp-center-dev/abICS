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
from abics.sampling.mc import verylargeint, MCAlgorithm, write_obs_header
from abics.sampling.mc_mpi import ParallelMC
from abics.util import pickle_dump, pickle_load, numpy_save, numpy_load


class RXParams:
    """Parameter set for replica exchange Monte Carlo

    Attributes
    ----------
    nreplicas : int
        The number of replicas
    nprocs_per_replica : int
        The number of processes which a replica uses
    kTs: np.ndarray[float]
        The temperature list
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
        self.nreplicas = 1
        self.nprocs_per_replica = 1
        self.kTs = np.zeros(0)
        self.nsteps = 0
        self.RXtrial_frequency = 1
        self.sample_frequency = 1
        self.print_frequency = 1
        self.throw_out = 0.5
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
        if "kTs" in d:
            params.kTs = np.array(d["kTs"])
            if params.kTs.size != params.nreplicas:
                raise ValueError(
                    "The number of temperatures must be equal to the number of replicas."
                )
        else:
            kTstart = d["kTstart"]
            kTend = d["kTend"]
            if d.get("linspace_in_beta", False):
                params.kTs = 1.0 / np.linspace(1.0 / kTstart, 1.0 / kTend, params.nreplicas)
            else:
                params.kTs = np.linspace(kTstart, kTend, params.nreplicas)
        params.nsteps = d["nsteps"]
        params.RXtrial_frequency = d.get("RXtrial_frequency", 1)
        params.sample_frequency = d.get("sample_frequency", 1)
        params.print_frequency = d.get("print_frequency", 1)
        params.reload = d.get("reload", False)
        params.seed = d.get("seed", 0)
        params.throw_out = d.get("throw_out", 0.5)
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


class TemperatureRX_MPI(ParallelMC):
    def __init__(
        self,
        comm,
        MCalgo: type[MCAlgorithm],
        model: Model,
        configs,
        kTs: list[float],
        write_node: bool=True,
        T2E: float=1.0,
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
        super().__init__(comm, MCalgo, model, configs, kTs, write_node=write_node, T2E=T2E)
        self.mycalc.kT = self.kTs[self.rank]
        self.mycalc.config = configs[self.rank]
        self.betas = 1.0 / np.array(self.kTs)
        self.rank_to_T = np.arange(0, self.procs, 1, dtype=np.int64)
        self.float_buffer = np.array(0.0, dtype=np.float64)
        self.int_buffer = np.array(0, dtype=np.int64)
        self.obs_save = []
        self.Trank_hist = []
        self.kT_hist = []
        self.Lreload = False
        self.ntrials = np.zeros(self.procs, dtype=np.int64)
        self.naccepted = np.zeros(self.procs, dtype=np.int64)
        if not (self.procs == self.nreplicas == len(self.kTs)):
            if self.rank == 0:
                print(
                    "ERROR: You have to set the number of replicas equal to the "
                    + "number of temperatures equal to the number of processes"
                )
            sys.exit(1)

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

    def Xtrial(self, XCscheme):
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
            self.comm.Send(np.array([self.mycalc.energy]), dest=exchange_rank, tag=1)
            self.comm.Recv([self.int_buffer, 1, MPI.INT], source=exchange_rank, tag=2)
            self.rank_to_T[self.rank] = self.int_buffer
        self.comm.Allgather(self.rank_to_T[self.rank], self.rank_to_T)
        self.mycalc.kT = self.kTs[self.rank_to_T[self.rank]]
        return

    def run(
        self,
        nsteps: int,
        RXtrial_frequency: int,
        sample_frequency: int = verylargeint,
        print_frequency: int = verylargeint,
        nsubsteps_in_step: int = 1,
        throw_out: int | float = 0.5,
        observer: ObserverBase = ObserverBase(),
        subdirs: bool = True,
        save_obs: bool = True,
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
        nsubsteps_in_step: int
            The number of Monte Carlo substeps in one MC step.
        throw_out: int | float
            The number (int) or ratio (float) of measurements to be dropped out as thermalization
        observer: observer object
        subdirs: boolean
            If true, working directory for this rank is made
        save_obs: boolean

        Returns
        -------
        obs_list: list
            Observation list
        """

        self.obsnames = observer.names

        if subdirs:
            try:
                os.mkdir(str(self.rank))
            except FileExistsError:
                pass
            os.chdir(str(self.rank))
        self.accept_count = 0
        if not self.Lreload:
            # self.mycalc.config.shuffle()
            self.mycalc.energy = self.mycalc.model.energy(self.mycalc.config)
        with open(os.devnull, "w") as f:
            test_observe = observer.observe(self.mycalc, f, lprint=False)
            obs_len = len(test_observe)
            obs = np.zeros([len(self.kTs), obs_len])
        # with open("obs.dat", "a") as output:
        #    obs_step = observer.observe(self.mycalc, output, True)

        nsample = 0
        XCscheme = 0
        with open("obs.dat", "a") as output:
            write_obs_header(output, self.mycalc, observer)
            ntrials = self.mycalc.ntrials
            naccepted = self.mycalc.naccepted
            for i in range(1, nsteps + 1):
                self.mycalc.MCstep(nsubsteps_in_step)
                if i % RXtrial_frequency == 0:
                    iT = self.rank_to_T[self.rank]
                    self.ntrials[iT] += self.mycalc.ntrials - ntrials
                    self.naccepted[iT] += self.mycalc.naccepted - naccepted
                    ntrials = self.mycalc.ntrials
                    naccepted = self.mycalc.naccepted

                    self.Xtrial(XCscheme)
                    XCscheme = (XCscheme + 1) % 2
                if i % sample_frequency == 0:
                    to_print = i % print_frequency == 0 and self.write_node
                    obs_step = observer.observe(self.mycalc, output, to_print)
                    obs[self.rank_to_T[self.rank]] += obs_step
                    if save_obs:
                        self.obs_save.append(obs_step)
                        self.Trank_hist.append(self.rank_to_T[self.rank])
                        self.kT_hist.append(self.mycalc.kT)
                    if self.write_node:
                        self.save(
                            save_obs=save_obs,
                            subdirs=subdirs,
                        )
                    nsample += 1
            iT = self.rank_to_T[self.rank]
            self.ntrials[iT] += self.mycalc.ntrials - ntrials
            self.naccepted[iT] += self.mycalc.naccepted - naccepted

        with open("acceptance_ratio.dat", "w") as f:
            for T, acc, trial in zip(self.kTs, self.naccepted, self.ntrials):
                if trial > 0:
                    f.write(f"{self.E2T*T} {acc/trial}\n")
                else:
                    f.write(f"{self.E2T*T} {acc}/{trial}\n")

        if nsample != 0:
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
                self.postproc(throw_out)
            return obs_list

        if subdirs:
            os.chdir("../")

    def save(self, save_obs: bool, subdirs: bool):
        self.comm.Barrier()
        # save information for restart
        pickle_dump(self.mycalc.config, "calc.pickle")
        rand_state = rand.get_state()
        pickle_dump(rand_state, "rand_state.pickle")
        if save_obs:
            obs_save_, Trank_hist_, kT_hist_ = self.__merge_obs()
            numpy_save(obs_save_, "obs_save.npy")
            numpy_save(Trank_hist_, "Trank_hist.npy")
            numpy_save(kT_hist_, "kT_hist.npy")
        if self.rank == 0:
            if subdirs:
                os.chdir("../")
            pickle_dump(self.rank_to_T, "rank_to_T.pickle")
            numpy_save(self.kTs, "kTs.npy")
            if subdirs:
                os.chdir(str(self.rank))

    def postproc(self, throw_out: int | float):
        assert throw_out >= 0
        obs_save, Trank_hist, kT_hist = self.__merge_obs()
        kTs = self.kTs
        comm = self.comm
        obsnames = self.obsnames
        postproc(obs_save, Trank_hist, kT_hist, kTs, comm, obsnames, throw_out, E2T=self.E2T)

    def __merge_obs(self):
        if hasattr(self, "obs_save0"):
            obs_save_ = np.concatenate((self.obs_save0, np.array(self.obs_save)))
            Trank_hist_ = np.concatenate((self.Trank_hist0, np.array(self.Trank_hist)))
            kT_hist_ = np.concatenate((self.kT_hist0, np.array(self.kT_hist)))
        else:
            obs_save_ = np.array(self.obs_save)
            Trank_hist_ = np.array(self.Trank_hist)
            kT_hist_ = np.array(self.kT_hist)
        return obs_save_, Trank_hist_, kT_hist_

def jackknife(X: np.ndarray) -> np.ndarray:
    nsamples = X.shape[0]
    S = X.sum(axis=0)
    X_jk = (S - X) / (nsamples - 1)
    return X_jk


def postproc(obs_save, Trank_hist, kT_hist, kTs, comm,
             obsnames, throw_out: int | float,
             E2T: float = 1.0,
             ):
    assert throw_out >= 0
    rank = comm.Get_rank()
    nT = comm.Get_size()
    nsteps, nobs = obs_save.shape
    # nT = rank_to_T.size
    if isinstance(throw_out, float):
        throw_out = int(nsteps * throw_out)

    # gather observables with same temperature (Trank == myrank)
    recv_obss = []
    start = throw_out
    buffer_size = 1000
    while start < nsteps:
        end = start + buffer_size
        if end > nsteps:
            end = nsteps
        obs = [[] for _ in range(nT)]
        for istep in range(start, end):
            iT = Trank_hist[istep]
            obs[iT].append(obs_save[istep, :])
        start = end
        send_obs = []
        for rk in range(comm.size):
            if len(obs[rk]) > 0:
                send = np.array(obs[rk])
            else:
                send = np.zeros((0, nobs))
            send_obs.append(send)
        recv_obs = comm.alltoall(send_obs)
        for ro in recv_obs:
            recv_obss.append(ro)
    X = np.concatenate(recv_obss)  # nsamples X nobs
    nsamples = X.shape[0]

    # jackknife method
    X2 = X**2
    X_mean = X.mean(axis=0)
    X_err = np.sqrt(X.var(axis=0, ddof=1) / (nsamples - 1))
    X_jk = jackknife(X)
    X2_mean = X2.mean(axis=0)
    X2_err = np.sqrt(X2.var(axis=0, ddof=1) / (nsamples - 1))
    X2_jk = jackknife(X2)
    F = X2.mean(axis=0) - X.mean(axis=0) ** 2  # F stands for Fluctuation
    F_jk = X2_jk - X_jk**2
    F_mean = nsamples * F - F_jk.mean(axis=0) * (nsamples - 1)
    F_err = np.sqrt((nsamples - 1) * F_jk.var(axis=0, ddof=0))

    # estimate free energy
    target_T = -1.0
    if kTs[0] <= kTs[-1]:
        if rank > 0:
            target_T = kTs[rank - 1]
    else:
        if rank < comm.size - 1:
            target_T = kTs[rank + 1]
    if target_T >= 0.0:
        dbeta = 1.0 / target_T - 1.0 / kTs[rank]
        energies = X[:, 0]
        emin = energies.min()
        Exp = np.exp(-dbeta * (energies - emin))
        Exp_jk = jackknife(Exp)
        Logz_jk = np.log(Exp_jk)
        Logz_mean = nsamples * np.log(Exp.mean()) - Logz_jk.mean() * (nsamples - 1) - dbeta * emin
        Logz_err = np.sqrt((nsamples - 1) * Logz_jk.var(ddof=0))
    else:
        Logz_mean = 0.0
        Logz_err = 0.0

    obs = np.array([X_mean, X_err, X2_mean, X2_err, F_mean, F_err])
    obs_all = np.zeros([comm.size, *obs.shape])  # nT X nobs X ntype
    comm.Allgather(obs, obs_all)
    dlogz = comm.allgather([Logz_mean, Logz_err])  # nT X 2

    if rank == 0:
        ntype = obs.shape[0]
        for iobs, oname in enumerate(obsnames):
            with open(f"{oname}.dat", "w") as f:
                f.write("# $1: temperature\n")
                f.write(f"# $2: <{oname}>\n")
                f.write(f"# $3: ERROR of <{oname}>\n")
                f.write(f"# $4: <{oname}^2>\n")
                f.write(f"# $5: ERROR of <{oname}^2>\n")
                f.write(f"# $6: <{oname}^2> - <{oname}>^2\n")
                f.write(f"# $7: ERROR of <{oname}^2> - <{oname}>^2\n")
                for iT in range(nT):
                    f.write(f"{E2T*kTs[iT]}")
                    for itype in range(ntype):
                        f.write(f" {obs_all[iT, itype, iobs]}")
                    f.write("\n")

        with open("logZ.dat", "w") as f:
            f.write("# $1: temperature\n")
            f.write("# $2: logZ\n")
            f.write("# $3: ERROR of log(Z)\n")
            f.write("# $4: log(Z/Z')\n")
            f.write("# $5: ERROR of log(Z/Z')\n")
            F = 0.0
            dF = 0.0
            if kTs[0] <= kTs[-1]:
                f.write(f"{E2T*kTs[-1]} {F} {dF} {0.0} {0.0}\n")
                for iT in np.arange(1, nT)[::-1]:
                    dlz, dlz_e = dlogz[iT]
                    F += dlz
                    dF += dlz_e
                    f.write(f"{E2T*kTs[iT-1]} {F} {dF} {dlz} {dlz_e}\n")
            else:
                f.write(f"{E2T*kTs[0]} {F} {dF} {0.0} {0.0}\n")
                for iT in np.arange(0, nT - 1):
                    dlz, dlz_e = dlogz[iT]
                    F += dlz
                    dF += dlz_e
                    f.write(f"{E2T*kTs[iT+1]} {F} {dF} {dlz} {dlz_e}\n")
    comm.Barrier()
