import os.path

import numpy as np


def read_coords(v):
    if isinstance(v, str):
        return np.array([[float(x) for x in line.split()] for line in v.splitlines()])
    else:
        return np.array(v)


def expand_path(path, basedir):
    path = os.path.expanduser(path)
    path = os.path.expandvars(path)
    if not path.startswith('/'):
        path = os.path.join(basedir, path)
    return path
