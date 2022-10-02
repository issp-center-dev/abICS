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

import sys

import numpy as np
import numpy.random as rand


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
