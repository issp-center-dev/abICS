from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
import numpy.random

from abics.mc import Model, ObserverBase, MCAlgorithm


class Configuration:
    Q: int
    spins: NDArray[np.int64]
    energy: float
    N: int
    N_zero: int

    def __init__(self, Q: int, Ls) -> None:
        self.Q = Q
        self.spins = np.zeros(Ls, dtype=np.int64)
        self.N = self.spins.size
        self.shuffle()

    def shuffle(self) -> None:
        self.spins = numpy.random.randint(self.Q, size=self.spins.shape)
        self.calc_energy()
        self.calc_mag()

    def calc_energy(self):
        self.energy = 0.0
        for index in np.ndindex(self.spins.shape):
            spin: int = self.spins[tuple(index)]
            npara = np.count_nonzero(self.neighbor_spins(index) == spin)
            self.energy -= npara
        self.energy *= 0.5  # due to double-counting

    def calc_mag(self):
        self.N_zero = self.N - np.count_nonzero(self.spins)

    def flip_spin(self, newspin, index, denergy=None) -> None:
        if denergy is None:
            denergy = self.diff_energy(newspin, index)
        self.energy += denergy
        if self.spins[index] == 0:
            self.N_zero -= 1
        self.spins[index] = newspin
        if self.spins[index] == 0:
            self.N_zero += 1

    def diff_energy(self, newspin, index) -> float:
        oldspin = self.spins[tuple(index)]
        neighbors = self.neighbor_spins(tuple(index))
        npara_old = np.count_nonzero(neighbors == oldspin)
        npara_new = np.count_nonzero(neighbors == newspin)
        return npara_old - npara_new

    def neighbor_spins(self, index) -> NDArray[np.int64]:
        ret = np.zeros(2*self.spins.ndim, dtype=np.int64)
        index_neighbor = np.array(index)
        for d in range(self.spins.ndim):
            index_neighbor[d] = (index[d] + 1) % self.spins.shape[d]
            ret[2*d] = self.spins[tuple(index_neighbor)]
            index_neighbor[d] = (index[d] - 1 + self.spins.shape[d]) % self.spins.shape[d]
            ret[2*d+1] = self.spins[tuple(index_neighbor)]
            index_neighbor[d] = index[d]
        return ret

    @property
    def mag(self) -> float:
        return self.N_zero / self.N - 1.0 / self.Q


@dataclass
class DConfig:
    newspin: int
    index: tuple[int, ...]


class Observer(ObserverBase):
    def logfunc(self, calc_state: MCAlgorithm):
        config: Configuration = calc_state.config
        spins = config.spins
        N = spins.size
        N_nonzero = np.count_nonzero(spins)
        N_zero = N - N_nonzero
        mag = N_zero / N - 1.0 / config.Q
        return (calc_state.energy, mag)


class Potts(Model):
    def __init__(self):
        ...

    def energy(self, config: Configuration) -> float:
        return config.energy

    def trialstep(self, config: Configuration, energy: float) -> tuple[DConfig, float]:
        index: tuple[int, ...] = tuple(np.random.randint(config.spins.shape))
        newspin = (config.spins[index] + np.random.randint(1, config.Q)) % config.Q
        denergy = config.diff_energy(newspin, index)
        return DConfig(newspin, index), denergy

    def newconfig(self, config: Configuration, dconfig: DConfig) -> Configuration:
        config.flip_spin(dconfig.newspin, dconfig.index)
        return config
