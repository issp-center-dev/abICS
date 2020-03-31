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

import os.path
import pickle

import numpy as np

from .exception import InputError


def read_vector(v, *, dtype=np.float):
    return read_tensor(v, rank=1, dtype=dtype)


def read_matrix(v, *, dtype=np.float):
    return read_tensor(v, rank=2, dtype=dtype)


def read_tensor(v, *, rank=2, dtype=np.float):
    """
    Read tensor

    Parameters
    ----------
    v: str or list or np.ndarray
        tensor

    dtype: type
        type of elements, default: np.float

    Returns
    -------
    v: numpy array
        tensor
    """

    def parse(x, dtype):
        if dtype == bool:
            if isinstance(x, str):
                xx = x.lower()
                if xx == "t" or xx == "true":
                    return True
                if xx == "f" or xx == "false":
                    return False
                n = int(x)
                return bool(n)
        else:
            return dtype(x)

    if isinstance(v, str):
        if rank == 1:
            ret = []
            lines = list(filter(lambda l: bool(l.strip()), v.splitlines()))
            n = len(lines)
            if n > 1:
                ret = [parse(x, dtype) for x in lines]
            else:
                ret = [parse(x, dtype) for x in lines[0].split()]
            return ret
        elif rank == 2:
            ret = []
            m0 = -1
            for line in v.splitlines():
                try:
                    row = [parse(x, dtype) for x in line.strip().split()]
                except ValueError as e:
                    raise InputError(str(e))

                if not row:
                    continue
                m = len(row)
                if m0 == -1:
                    m0 = m
                if m != m0:
                    raise InputError("Dimension mismatch in {}".format(v))
                ret.append(row)
            if not ret:
                return np.zeros((0, 0), dtype=dtype)
            return np.array(ret, dtype=dtype)
        else:
            raise InputError("read_tensor(rank>2) requires list-type argument")
    elif isinstance(v, list):
        if rank == 1:
            return np.array(v, dtype=dtype)
        m0 = -1
        ret = []
        for vv in v:
            child = read_tensor(vv, rank=rank-1, dtype=dtype)
            m = child.shape[0]
            if m0 < 0:
                m0 = m
            if m != m0:
                raise InputError("Dimension mismatch in {}".format(v))
            ret.append(child)
        return np.array(ret, dtype=dtype)
    else:
        return np.array(v, dtype=dtype)


def expand_path(path, basedir):
    """
    Expand path (absolute path)

    Parameters
    ----------
    path: str
        path

    basedir: str
        path to base directory

    Returns
    -------
    path: str
        path combined basedir with path (absolute path)
    """
    path = os.path.expanduser(path)
    path = os.path.expandvars(path)
    if not path.startswith('/'):
        path = os.path.join(basedir, path)
    return path


def pickle_dump(data, filename):
    with open(filename, 'wb') as f:
        pickle.dump(data, f)


def pickle_load(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

def numpy_save(data, filename, allow_pickle=False):
    with open(filename, 'wb') as f:
        np.save(f, data, allow_pickle)


def numpy_load(filename):
    with open(filename, 'rb') as f:
        return np.load(f)
