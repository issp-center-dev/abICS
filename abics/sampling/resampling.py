# original file is distributed from py2dmat: https://github.com/issp-center-dev/2DMAT
# under the GNU GPL v3 (see the following)

# 2DMAT -- Data-analysis software of quantum beam diffraction experiments for 2D material structure
# Copyright (C) 2020- The University of Tokyo
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

import typing
from typing import Union, List, Iterable

import abc

import collections
import itertools
import numpy as np


class Resampler(abc.ABC):
    @abc.abstractmethod
    def reset(self, weights: Iterable):
        ...

    @abc.abstractmethod
    def sample(self, size=None) -> Union[int, np.ndarray]:
        ...


class BinarySearch(Resampler):
    weights_accumulated: List[float]
    wmax: float

    def __init__(self, weights: Iterable):
        self.reset(weights)

    def reset(self, weights: Iterable):
        self.weights_accumulated = list(itertools.accumulate(weights))
        self.wmax = self.weights_accumulated[-1]

    @typing.overload
    def sample(self) -> int:
        ...

    @typing.overload
    def sample(self, size) -> np.ndarray:
        ...

    def sample(self, size=None) -> Union[int, np.ndarray]:
        if size is None:
            return self._sample(self.wmax * np.random.rand())
        else:
            return np.array([self._sample(r) for r in self.wmax * np.random.rand(size)])

    def _sample(self, r: float) -> int:
        return typing.cast(int, np.searchsorted(self.weights_accumulated, r))


class WalkerTable(Resampler):
    N: int
    itable: np.ndarray
    ptable: np.ndarray

    def __init__(self, weights: Iterable):
        self.reset(weights)

    def reset(self, weights: Iterable):
        self.ptable = np.array(weights).astype(np.float64).flatten()
        self.N = len(self.ptable)
        self.itable = np.full(self.N, -1)
        mean = self.ptable.mean()
        self.ptable -= mean
        shorter = collections.deque([i for i, p in enumerate(self.ptable) if p < 0.0])
        longer = collections.deque([i for i, p in enumerate(self.ptable) if p >= 0.0])

        while len(longer) > 0 and len(shorter) > 0:
            ilong = longer[0]
            ishort = shorter.popleft()
            self.itable[ishort] = ilong
            self.ptable[ilong] += self.ptable[ishort]
            if self.ptable[ilong] <= 0.0:
                longer.popleft()
                shorter.append(ilong)
        self.ptable += mean
        self.ptable *= 1.0 / mean

    @typing.overload
    def sample(self) -> int:
        ...

    @typing.overload
    def sample(self, size) -> np.ndarray:
        ...

    def sample(self, size=None) -> Union[int, np.ndarray]:
        if size is None:
            r = np.random.rand() * self.N
            return self._sample(r)
        else:
            r = np.random.rand(size) * self.N
            i = np.floor(r).astype(np.int64)
            p = r - i
            ret = np.where(p < self.ptable[i], i, self.itable[i])
            return ret

    def _sample(self, r: float) -> int:
        i = int(np.floor(r))
        p = r - i
        if p < self.ptable[i]:
            return i
        else:
            return self.itable[i]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="resampling test")
    parser.add_argument(
        "-s", "--seed", type=int, default=12345, help="random number seed"
    )
    parser.add_argument("-m", "--method", default="walker", help="method to resample")
    parser.add_argument("-N", type=int, default=100000, help="Number of samples")

    args = parser.parse_args()

    # rs = np.random.RandomState(args.seed)

    ps = [0.0, 1.0, 2.0, 3.0]
    # ps = rs.rand(5)
    S = np.sum(ps)
    if args.method == "walker":
        resampler = WalkerTable(ps)
    else:
        resampler = BinarySearch(ps)
    samples = resampler.sample(args.N)

    print("#i result exact diff")
    for i, p in enumerate(ps):
        r = np.count_nonzero(samples == i) / args.N
        print(f"{i} {r} {p/S} {r - p/S}")
