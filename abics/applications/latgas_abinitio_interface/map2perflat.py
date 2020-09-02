#from mpi4py import MPI
import os
import numpy as np
from pymatgen import Structure
from abics.applications.latgas_abinitio_interface import aenet, naive_matcher

mapper = naive_matcher.naive_mapping

def map2perflat(perf_st, st, vac_spaceholder={}):
    perf_st = perf_st.copy()
    N = perf_st.num_sites
    seldyn = np.array(perf_st.site_properties.get("seldyn", np.ones((N, 3))), dtype=np.float)
    perf_st.replace_species(vac_spaceholder)        
    mapping = mapper(st, perf_st)
    for i,j in enumerate(mapping):
        perf_st.replace(j, st[i].species_string, properties={"seldyn":seldyn[j]})
    perf_st.sort(key=lambda site: site.species_string)
    return perf_st

if __name__ == '__main__':
    perf_st = Structure.from_file("../../POSCAR")
    st_in = Structure.from_file("./POSCAR")
    st = map2perflat(perf_st, st_in)
    st.remove_species(["N", "Li"])
    st.to("POSCAR", "POSCAR.in")
