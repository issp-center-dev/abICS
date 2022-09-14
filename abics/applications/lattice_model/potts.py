from typing import Any, List, Tuple

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
import numpy.random

from abics.mc import Model


class Configuration:
    Q: int
    spins: NDArray[np.int64]
    energy: float

    def __init__(self, Q: int, Ls) -> None:
        self.Q = Q
        self.spins = np.zeros(Ls, dtype=np.int64)
        self.shuffle()

    def shuffle(self) -> None:
        self.spins = numpy.random.randint(self.Q, size=self.spins.shape)
        self.calc_energy()

    def calc_energy(self):
        self.energy = 0.0
        for index in np.ndindex(self.spins.shape):
            spin: int = self.spins[tuple(index)]
            npara = np.count_nonzero(self.neighbor_spins(index) == spin)
            self.energy -= npara

    def flip_spin(self, newspin, index, denergy=None) -> None:
        if denergy is None:
            denergy = self.diff_energy(newspin, index)
        self.energy += denergy
        self.spins[index] = newspin

    def diff_energy(self, newspin, index) -> float:
        oldspin = self.spins[tuple(index)]
        neighbors = self.neighbor_spins(index)
        npara_old = np.count_nonzero(neighbors == oldspin)
        npara_new = np.count_nonzero(neighbors == newspin)
        return npara_old - npara_new

    def neighbor_spins(self, index) -> NDArray[np.int64]:
        ret = np.zeros(self.spins.ndim, dtype=np.int64)
        index_neighbor = np.array(index)
        for d in range(self.spins.ndim):
            index_neighbor[d] = (index_neighbor[d] + 1) % self.spins.shape[d]
            ret[d] = self.spins[tuple(index_neighbor)]
            index_neighbor[d] = index[d]
        return ret


@dataclass
class DConfig:
    newspin: int
    index: NDArray[np.int64]


class Potts(Model):
    def __init__(self):
        ...

    def energy(self, config: Configuration) -> float:
        return config.energy

    def trialstep(self, config: Configuration, energy: float) -> Tuple[DConfig, float]:
        index = np.random.randint(config.spins.shape)
        newspin = (config.spins[index] + np.random.randint(1, config.Q)) % config.Q
        denergy = config.diff_energy(newspin, index)
        return DConfig(newspin, index), denergy

    def newconfig(self, config: Configuration, dconfig: DConfig):
        config.flip_spin(dconfig.newspin, dconfig.index)
        return config
