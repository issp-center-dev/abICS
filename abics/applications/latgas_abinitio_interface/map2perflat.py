from __future__ import annotations

import numpy as np
from pymatgen.core import Structure
from abics.applications.latgas_abinitio_interface.naive_matcher import naive_mapping


def map2perflat(
    perf_st: Structure, st: Structure, vac_spaceholder: dict = {}
) -> Structure:
    """

    Arguments
    =========
    perf_st: Structure
        reference structure
    st: Structure
        input structure (configuration)
    vac_spaceholder: dict[str, str]
        - key
            - specy 'A' maybe removed
        - value
            - virtual specy implying vacancy of 'A'

    """
    perf_st = perf_st.copy()
    N = perf_st.num_sites
    seldyn = np.array(
        perf_st.site_properties.get("seldyn", np.ones((N, 3))), dtype=np.float64
    )
    perf_st.replace_species(vac_spaceholder)
    mapping = naive_mapping(st, perf_st)
    for i, j in enumerate(mapping):
        perf_st.replace(j, st[i].species_string, properties={"seldyn": seldyn[j]})
    perf_st.sort(key=lambda site: site.species_string)
    return perf_st


"""
Testing
"""
if __name__ == "__main__":
    perf_st = Structure.from_file("../../POSCAR")
    st_in = Structure.from_file("./POSCAR")
    st = map2perflat(perf_st, st_in)
    st.remove_species(["N", "Li"])
    st.to(fmt="POSCAR", filename="POSCAR.in")
