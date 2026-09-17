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

from typing import Any, Optional

from abc import ABCMeta, abstractmethod

# import SFMT_cython.sfmt_random as sfmt_random
# from multiprocessing import Process, Queue, Pool, TimeoutError
import os
import sys
import numpy as np
import numpy.random as rand

from abics import __version__
from abics.observer import ObserverBase
from abics.model import Model
from abics.sampling.mc import MCAlgorithm

import logging

logger = logging.getLogger("main")

verylargeint = sys.maxsize


def binning(x, nlevels: int):
    """

    Parameters
    ----------
    x: list
        Coordinates
    nlevels: int
        Number to determine the number of measurements

    Returns
    -------
    error_estimate: list
        error estimation
    """
    error_estimate = []
    x = np.array(x, dtype=np.float64)
    # assert 2**nlevels*10 < len(x)
    throwout = len(x) % (2**nlevels)
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





class WangLandauMonteCarlo(MCAlgorithm):
    def __init__(
        self,
        model: Model,
        finit: float,
        emin: float,
        emax: float,
        n_interval: int,
        config,
        ffactor: float = 0.9,
        flatness: float = 0.8,
        check_interval: int = 1000,
        ffinal: float = 0.0,
    ):
        super().__init__()

        if emin >= emax:
            raise ValueError("emin must be smaller than emax")
        if finit <= 0.0:
            raise ValueError("finit must be positive")
        if ffinal < 0.0:
            raise ValueError("ffinal must be non-negative")
        if ffinal > 0.0 and ffinal >= finit:
            raise ValueError("ffinal must be smaller than finit")
        if not 0.0 < ffactor < 1.0:
            raise ValueError("ffactor must be between zero and one")
        if not 0.0 < flatness <= 1.0:
            raise ValueError("flatness must be in (0, 1]")
        if check_interval < 1:
            raise ValueError("check_interval must be positive")
        if n_interval < 2:
            raise ValueError("n_interval must be at least 2")

        self.model = model
        self.config = config
        self.emin = float(emin)
        self.emax = float(emax)
        self.n_interval = int(n_interval)
        self.bin_width = (self.emax - self.emin) / self.n_interval
        self.energies = self.emin + (
            np.arange(self.n_interval) + 0.5
        ) * self.bin_width

        self.finit = float(finit)
        self.fcurrent = float(finit)
        self.ffinal = float(ffinal)
        self.ffactor = float(ffactor)
        self.flatness = float(flatness)
        self.check_interval = int(check_interval)
        self._steps_since_check = 0

        self.hist = np.zeros(self.n_interval, dtype=np.int64)
        self.entropy = np.zeros(self.n_interval, dtype=np.float64)
        # bins ever visited across the whole run; unlike self.hist, never reset
        self._ever_visited = np.zeros(self.n_interval, dtype=bool)
        self.obs_save = []

        logger.info("Wang Landau: finit = %s", self.finit)
        logger.info("Wang Landau: emin = %s", self.emin)
        logger.info("Wang Landau: emax = %s", self.emax)
        logger.info("Wang Landau: n_interval = %s", self.n_interval)

    def initialize(self) -> None:
        self.energy = self.model.energy(self.config)
        if self._energy_index(self.energy) is None:
            raise ValueError(
                f"Initial energy {self.energy} is outside "
                f"[{self.emin}, {self.emax})"
            )

    def finished(self) -> bool:
        return self.ffinal > 0.0 and self.fcurrent <= self.ffinal

    def _energy_index(self, energy: float) -> Optional[int]:
        if energy < self.emin or energy >= self.emax:
            return None
        index = int((energy - self.emin) / self.bin_width)
        return min(index, self.n_interval - 1)

    def _check_flatness(self) -> None:
        # bins never visited (e.g. a true gap in the DOS) are excluded so they
        # cannot block convergence forever
        if not np.any(self._ever_visited):
            return
        if np.any(self.hist[self._ever_visited] == 0):
            return

        visited_hist = self.hist[self._ever_visited]
        histogram_flatness = visited_hist.min() / visited_hist.mean()
        if histogram_flatness >= self.flatness:
            logger.info(
                "Histogram is sufficiently flat: %.3f",
                histogram_flatness,
            )
            self.fcurrent *= self.ffactor
            self.hist.fill(0)

    def MCstep(self, nsubsteps_in_step: int = 1):
        for _ in range(nsubsteps_in_step):
            current_index = self._energy_index(self.energy)
            if current_index is None:
                raise ValueError(
                    f"Current energy {self.energy} is outside "
                    f"[{self.emin}, {self.emax})"
                )

            dconfig, dE = self.model.trialstep(self.config, self.energy)
            proposed_energy = self.energy + dE
            proposed_index = self._energy_index(proposed_energy)

            accepted = False
            if proposed_index is not None:
                log_probability = (
                    self.entropy[current_index]
                    - self.entropy[proposed_index]
                )
                accepted = np.log(rand.random()) < min(
                    0.0, log_probability
                )

            if accepted:
                self.config = self.model.newconfig(self.config, dconfig)
                self.energy = proposed_energy
                visited_index = proposed_index
                self.naccepted += 1
            else:
                visited_index = current_index

            self.hist[visited_index] += 1
            self._ever_visited[visited_index] = True
            self.entropy[visited_index] += self.fcurrent
            self.ntrials += 1
            self._steps_since_check += 1

            if self._steps_since_check >= self.check_interval:
                self._check_flatness()
                self._steps_since_check = 0

            logger.debug(
                "MCstep: %s", "accepted" if accepted else "rejected"
            )

    def save_entropy(self, filename: str = "entropy.dat"):
        np.savetxt(
            filename,
            np.column_stack((self.energies, self.entropy)),
            header="Energy Entropy",
        )

    def save_dos(self, filename: str = "dos.dat"):
        dos = np.exp(self.entropy - np.max(self.entropy))
        np.savetxt(
            filename,
            np.column_stack((self.energies, dos)),
            header="Energy DOS",
        )

    def run(
        self,
        nsteps: int,
        sample_frequency: int = verylargeint,
        print_frequency: int = verylargeint,
        nsubsteps_in_step: int = 1,
        observer: ObserverBase = ObserverBase(),
        save_obs: bool = False,
    ):
        logger.info("WangLandauMonteCarlo: run")
        result = super().run(
            nsteps,
            sample_frequency,
            print_frequency,
            nsubsteps_in_step,
            observer,
            save_obs,
        )
        self.save_dos()
        self.save_entropy()
        return result

    def parameters(self):
        return [self.fcurrent]

    def param_names(self):
        return ["Modification factor"]


def write_obs_header(output, calc_state: MCAlgorithm, observer: ObserverBase):
    output.write("# $1: step\n")
    i = 2
    for name in calc_state.param_names():
        output.write(f"# ${i}: {name}\n")
        i+=1
    for name in observer.obs_names():
        output.write(f"# ${i}: {name}\n")
        i+=1


