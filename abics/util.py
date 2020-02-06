import os.path
import pickle

import numpy as np


def read_coords(v):
    """
    Read coordinates

    Parameters
    ----------
    v: str or numpy array
        coordinates information

    Returns
    -------
    v: numpy array
        coordinates information
    """
    if isinstance(v, str):
        return np.array([[float(x) for x in line.split()] for line in v.splitlines()])
    else:
        return np.array(v)


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
