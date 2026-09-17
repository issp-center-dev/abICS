from __future__ import annotations

import os

from mpi4py import MPI

import numpy as np
import numpy.random as rand

from abics.model import Model
from abics.observer import ObserverBase
from abics.sampling.mc import verylargeint
from abics.sampling.wl import WangLandauMonteCarlo
from abics.util import pickle_dump, pickle_load


class ParallelWL:
    """Run independent Wang-Landau walkers on MPI ranks.

    Each rank owns one configuration and samples the same energy window.
    Entropy and histogram arrays remain rank-local; this class runs
    independent walkers rather than combining their density-of-states data.
    """

    comm: MPI.Comm
    rank: int
    procs: int
    model: Model
    mycalc: WangLandauMonteCarlo
    energy_window: tuple[float, float]
    write_node: bool
    Lreload: bool

    def __init__(
        self,
        comm,
        MCalgo: type[WangLandauMonteCarlo],
        model: Model,
        configs,
        energy_window: tuple[float, float],
        n_interval: int,
        finit: float,
        ffactor: float = 0.9,
        flatness: float = 0.8,
        check_interval: int = 1000,
        write_node: bool = True,
    ):
        self.comm = comm
        self.rank = comm.Get_rank()
        self.procs = comm.Get_size()
        self.model = model
        self.write_node = write_node
        self.Lreload = False

        if len(configs) != self.procs:
            raise ValueError("configs must contain one configuration per MPI rank")
        if n_interval < 2:
            raise ValueError("n_interval must be at least 2")

        if len(energy_window) != 2:
            raise ValueError("energy_window must contain (emin, emax)")
        emin, emax = map(float, energy_window)
        if emin >= emax:
            raise ValueError("energy_window must have emin < emax")
        self.energy_window = (emin, emax)
        self.mycalc = MCalgo(
            model,
            finit,
            emin,
            emax,
            n_interval,
            configs[self.rank],
            ffactor=ffactor,
            flatness=flatness,
            check_interval=check_interval,
        )

    @property
    def entropy(self):
        return self.mycalc.entropy

    @property
    def histogram(self):
        return self.mycalc.hist

    def reload(self, subdirs: bool = True) -> None:
        """Reload the rank-local state saved by :meth:`run`."""
        state_dir = str(self.rank) if subdirs else "."
        state = pickle_load(os.path.join(state_dir, "wl_state.pickle"))
        self.mycalc.config = state["config"]
        self.mycalc.energy = state["energy"]
        self.mycalc.entropy = state["entropy"]
        self.mycalc.hist = state["hist"]
        self.mycalc._ever_visited = state["ever_visited"]
        self.mycalc.fcurrent = state["fcurrent"]
        self.mycalc.naccepted = state["naccepted"]
        self.mycalc.ntrials = state["ntrials"]
        self.mycalc._steps_since_check = state["steps_since_check"]
        rand.set_state(state["random_state"])
        self.Lreload = True

    def _save_state(self, filename: str) -> None:
        state = {
            "config": self.mycalc.config,
            "energy": self.mycalc.energy,
            "entropy": self.mycalc.entropy,
            "hist": self.mycalc.hist,
            "ever_visited": self.mycalc._ever_visited,
            "fcurrent": self.mycalc.fcurrent,
            "naccepted": self.mycalc.naccepted,
            "ntrials": self.mycalc.ntrials,
            "steps_since_check": self.mycalc._steps_since_check,
            "random_state": rand.get_state(),
        }
        pickle_dump(state, filename)

    def run(
        self,
        nsteps: int,
        sample_frequency: int = verylargeint,
        print_frequency: int = verylargeint,
        nsubsteps_in_step: int = 1,
        observer: ObserverBase = ObserverBase(),
        subdirs: bool = True,
    ):
        """Run one Wang-Landau walker per MPI rank.

        Returns an array of rank-local averaged observables when samples were
        collected. Entropy and DOS files are written by each rank locally.
        """
        if subdirs:
            os.makedirs(str(self.rank), exist_ok=True)
            os.chdir(str(self.rank))

        observables = self.mycalc.run(
            nsteps,
            sample_frequency=sample_frequency,
            print_frequency=print_frequency,
            nsubsteps_in_step=nsubsteps_in_step,
            observer=observer,
        )

        if self.write_node:
            pickle_dump(self.mycalc.config, "config.pickle")
            self._save_state("wl_state.pickle")

        if subdirs:
            os.chdir("../")

        self.comm.Barrier()
        if observables is None:
            return None

        observables = np.asarray(observables, dtype=float)
        obs_buffer = np.empty((self.procs, observables.size), dtype=float)
        self.comm.Allgather(observables, obs_buffer)
        return obs_buffer
