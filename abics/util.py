import os.path
import pickle

import numpy as np

from .exception import InputError

def read_matrix(v, *, dtype=np.float):
    """
    Read matrix

    Parameters
    ----------
    v: str or list of list or np.ndarray
        matrix

    dtype: type
        type of elements, default: np.float

    Returns
    -------
    v: numpy array
        matrix
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
            return np.zeros((0,0), dtype=dtype)
        return np.array(ret, dtype=dtype)
    elif isinstance(v, list):
        if not v:
            return np.zeros((0,0), dtype=dtype)
        m0 = -1
        for vv in v:
            if not isinstance(vv, list):
                raise InputError("{} is not list of list".format(v))
            m = len(vv)
            if m0 < 0:
                m0 = m
            if m != m0:
                raise InputError("Dimension mismatch in {}".format(v))
        return np.array(v, dtype=dtype)
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
