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

from typing import TextIO, Union, TYPE_CHECKING

import sys, os
import numpy as np

if TYPE_CHECKING:
    from mpi4py import MPI

from abics import __version__


class ObsInfo:
    nargs: int
    lengths: list[int]

    def __init__(self, *args):
        """

        Parameters
        ----------
        args: list
        """
        self.nargs = len(args)
        self.lengths = []
        for arg in args:
            # Inelegant way to make everything a 1D array
            arg = np.array([arg])
            arg = arg.ravel()
            self.lengths.append(len(arg))

    def decode(self, obs_array):
        """

        Parameters
        ----------
        obs_array: numpy array

        Returns
        -------
        obs: list

        """

        obs = []
        idx = 0
        for i in range(self.nargs):
            length = self.lengths[i]
            if length == 1:
                obs.append(obs_array[idx])
            else:
                obs.append(obs_array[idx : idx + length])
            idx += length
        return obs


# @profile
def obs_encode(*args):
    """make pure 1D data

    Parameters
    ----------
    args: list

    Returns
    -------
    obs_array: numpy array

    """
    # nargs = np.array([len(args)])
    # args_length_list = []
    obs_array = np.array([])
    for arg in args:
        # Inelegant way to make everything a 1D array
        arg = np.array([arg])
        arg = arg.ravel()
        obs_array = np.concatenate((obs_array, arg))
        # args_length_list.append(len(arg))
    # args_length_array = np.array(args_length_list)
    # args_info = np.concatenate((nargs, args_length_array))
    return obs_array


def obs_column_names(names, lengths=None) -> list[str]:
    """Expand names of array observables to one name per element, e.g. foo[0], foo[1]"""
    if lengths is None:
        lengths = [1] * len(names)
    cols = []
    for name, length in zip(names, lengths):
        if length == 1:
            cols.append(name)
        else:
            cols.extend(f"{name}[{i}]" for i in range(length))
    return cols


def write_obs_stats_file(oname, length, kTs, E2T, values, offset):
    """Write <oname>.dat of the temperature-dependent statistics.

    Parameters
    ----------
    oname: str
        observable name
    length: int
        number of elements of the observable
    kTs: array-like
        temperatures (kT)
    E2T: float
        conversion factor from kT to temperature
    values: callable
        values(iT, icol) -> sequence of 2*ntype values
        (mean, err, mean of square, err, fluctuation, err) for column icol
    offset: int
        index of the first column of this observable
    """
    with open(f"{oname}.dat", "w") as f:
        f.write("# $1: temperature\n")
        idx = 2
        for i in range(length):
            n = oname if length == 1 else f"{oname}[{i}]"
            f.write(f"# ${idx}: <{n}>\n")
            f.write(f"# ${idx+1}: ERROR of <{n}>\n")
            f.write(f"# ${idx+2}: <{n}^2>\n")
            f.write(f"# ${idx+3}: ERROR of <{n}^2>\n")
            f.write(f"# ${idx+4}: <{n}^2> - <{n}>^2\n")
            f.write(f"# ${idx+5}: ERROR of <{n}^2> - <{n}>^2\n")
            idx += 6
        for iT, kT in enumerate(kTs):
            f.write(f"{E2T*kT}")
            for i in range(length):
                for v in values(iT, offset + i):
                    f.write(f" {v}")
            f.write("\n")


def write_obs_npz(oname, length, offset, kTs, E2T, stats):
    """Write statistics of one observable into <oname>.npz.

    Parameters
    ----------
    oname: str
        observable name
    length: int
        number of elements of the observable
    offset: int
        index of the first column of this observable in `stats`
    stats: numpy array, shape (nT, 6, ncolumns)
        (mean, err, mean of square, err, fluctuation, err) for each column

    The file contains `temperature` (nT) and `mean`, `mean_err`,
    `second_moment`, `second_moment_err`, `fluctuation`, `fluctuation_err`.
    Each has shape (nT,) for a scalar observable and (nT, length) for an
    array observable.
    """
    stats = np.asarray(stats)
    keys = ["mean", "mean_err", "second_moment", "second_moment_err",
            "fluctuation", "fluctuation_err"]
    data = {}
    for i, k in enumerate(keys):
        v = stats[:, i, offset : offset + length]
        data[k] = v[:, 0] if length == 1 else v
    np.savez(f"{oname}.npz", temperature=E2T * np.asarray(kTs), **data)


class ObserverBase:
    comm: (None|MPI.Comm)
    Lreload: bool
    params: dict
    lprintcount: int
    names: list[str]

    def __init__(self, comm: (None|MPI.Comm) = None, Lreload: bool = False, params: dict = {}):
        self.comm = comm
        self.lprintcount = 0
        self.Lreload = Lreload
        self.params = params
        self.names = ["energy"]

    def obs_info(self, calc_state: "MCAlgorithm") -> ObsInfo:
        """

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm

        Returns
        -------
        args_info: numpy array

        """
        obs_log = self.logfunc(calc_state)
        if not isinstance(obs_log, tuple):
            obs_log = (obs_log,)
        obs_ND = []
        obs_save = self.savefunc(calc_state)
        if len(obs_save) > 0:
            if isinstance(obs_save, tuple):
                for obs in obs_save:
                    obs_ND.append(obs)
            else:
                obs_ND.append(obs_save)
        return ObsInfo(*obs_log, *obs_ND)

    def logfunc(self, calc_state: "MCAlgorithm") -> tuple[float, ...]:
        """returns values of observables

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm

        Returns
        -------
        calc_state.energy: tuple
            (Energy)
        """
        energy = calc_state.model.energy(calc_state.config)
        return (energy,)

    def savefunc(self, calc_state: "MCAlgorithm") -> tuple[float, ...]:
        """returns values of observables, which will be not printed in observe method.

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm

        Returns
        -------
        """
        return ()

    def writefile(self, calc_state: "MCAlgorithm") -> None:
        """

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm

        Returns
        -------

        """
        return None

    def observe(self, calc_state: "MCAlgorithm", outputfi: TextIO, lprint=True):
        """

        Parameters
        ----------
        calc_state: MCAlgorithm
            MonteCarlo algorithm
        outputfi: _io.TextIOWrapper
            TextIOWrapper for output
        lprint: boolean
            if true, log info is outputted to TextIOWrapper

        Returns
        -------
        obs_log: numpy array
            log information about observation
        """
        obs_log = self.logfunc(calc_state)
        if not isinstance(obs_log, tuple):
            obs_log = (obs_log,)
        obs_log = obs_encode(*obs_log)
        if lprint:
            line = f"{self.lprintcount}\t"
            for p in calc_state.parameters():
                line += f"{p}\t"
            for x in obs_log:
                line += f"{x}\t"
            line += "\n"
            outputfi.write(line)
            outputfi.flush()
            self.writefile(calc_state)
            self.lprintcount += 1
        obs_save = np.atleast_1d(self.savefunc(calc_state))
        if len(obs_save) > 0:
            obs_save = np.atleast_1d(obs_save)
            obs_save = obs_save.ravel()
            return np.concatenate((obs_log, obs_save))
        else:
            return obs_log

    def obs_names(self) -> list[str]:
        """

        Returns
        -------
        names: list
            names of observables
        """
        return self.names

    def obs_lengths(self) -> list[int]:
        """number of elements of each observable in `names`.

        Observers returning array-valued observables should override this.
        The lengths may be determined by the first call of `logfunc`.
        """
        return [1] * len(self.names)

    def obs_column_names(self) -> list[str]:
        """names of the columns written to obs.dat (array observables are expanded)"""
        return obs_column_names(self.obs_names(), self.obs_lengths())


class ObserverParams:
    type: str
    observer_class: (type | None)
    dict: dict

    def __init__(self):
        self.type = "default"
        self.observer_class = None
        self.dict = {}

    @classmethod
    def from_dict(cls, d: dict):
        """

        Parameters
        ----------
        d: dict
            Dictionary

        Returns
        -------
        oparams: ObserverParams
            self
        """
        params = cls()
        params.type = d.get("type", "default")
        if params.type == "import":
            sys.path.append(os.getcwd())
            from observer_module import Observer

            params.observer_class = Observer

        params.dict = {k: v for k, v in d.items() if k != "type"}
        return params

    @classmethod
    def from_toml(cls, f):
        """

        Parameters
        ----------
        f: str
            Name of input toml File

        Returns
        -------
        oparams : ObserverParams
            self
        """
        import toml

        d = toml.load(f)
        return cls.from_dict(d["observer"])


# For backward compatibility
if __version__ < "3":
    observer_base = ObserverBase
