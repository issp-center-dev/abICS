from typing import List

import numpy as np
from pymatgen.core import Structure


def naive_mapping(st0: Structure, st1: Structure) -> List[int]:
    """

    Arguments
    =========
    st0: Structure
    st1: Structure

    Returns
    =======
    mapping: List[int]
        mapping atoms from st0 to st1.
        st0[i] corresponds to st1[mapping[i]]
    """
    coords0 = st0.frac_coords
    coords1 = st1.frac_coords
    D = st1.lattice.get_all_distances(coords0, coords1)
    rows = D.min(axis=1).argsort()
    cols = D.argmin(axis=1)[rows]
    mapping = cols[rows.argsort()]
    return mapping


"""
Testing
"""
if __name__ == "__main__":
    st0 = Structure.from_file("MgAl2O4.vasp")
    st1 = Structure.from_file("0/structure.0.vasp")
    for i, j in enumerate(naive_mapping(st0, st1)):
        st0.replace(i, st1[j].species_string)
    st0.sort()
    st0.to("POSCAR", "testnaivemap.vasp")
