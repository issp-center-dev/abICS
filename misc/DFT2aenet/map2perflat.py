from mpi4py import MPI
import os
import numpy as np
from pymatgen import Structure
from abics.applications.latgas_abinitio_interface import aenet
import naive_matcher
mapper = naive_matcher.naive_mapping

def map2perflat(perf_st, st, vac_spaceholder={}):
    perf_st = perf_st.copy()
    perf_st.replace_species(vac_spaceholder)        
    mapping = mapper(st, perf_st)
    for i,j in enumerate(mapping):
        perf_st.replace(j, st[i].species_string)
    perf_st.sort()
    return perf_st


default_spaceholders = ["He", "Ne", "Ar", "Kr", "Xe", "Rn"]
vac_spaceholder = {"O":"N", "H":"Li"}

consider_only = ["Mg", "Al"]#["Sc", "H", "N"]
#["Zr", "O", "Li"]

comm = MPI.COMM_WORLD
myreplica = comm.Get_rank()

if myreplica == 0:
    os.makedirs("aenet_xsf", exist_ok = True)
comm.Barrier()

energies = np.loadtxt(os.path.join(str(myreplica),"energy_corr.dat"))[:,1] 
nstructure = energies.shape[0]
perf_st = Structure.from_file("../MgAl2O4.vasp")
perf_st.make_supercell([2,2,2])
perf_st.replace_species({"Al":"Mg"})
sublattice_mapping = {"Mg":["Mg", "Al"]}

structure_ids = list(range(0, 400, 20))
for i in range(nstructure):
    stmapped = []
    st = Structure.from_file(os.path.join(str(myreplica),"structure_rel.{}.vasp".format(structure_ids[i])))
    for sublattice in sublattice_mapping.keys():
        perf_st_tmp = perf_st.copy()
        perf_st_tmp.remove_species(
            [sp for sp in perf_st_tmp.symbol_set if sp != sublattice]
            )
        if sublattice in vac_spaceholder.keys():
            perf_st_tmp.replace_species({sublattice: vac_spaceholder[sublattice]})
        st_tmp = st.copy()
        st_tmp.remove_species(
            [sp for sp in st_tmp.symbol_set if sp not in sublattice_mapping[sublattice]]
            )
        #print(perf_st_tmp)
        #print(st_tmp)
        stmapped.append(map2perflat(perf_st_tmp, st_tmp))

    st_tmp = stmapped[0]
    for st1 in stmapped[1:]:
        for site in st1:
            st_tmp.append(site.species_string, site.frac_coords)
    st_tmp.remove_species(
            [sp for sp in st_tmp.symbol_set if sp not in consider_only]
            )
            
    #st.remove_species(remove_list)
    #perf_st.remove_species(remove_list)
    #stmapped = map2perflat(perf_st, st, vac_spaceholder)
    #sp = list(stmapped.symbol_set)
    #remove_after_mapping = [spec for spec in sp if spec not in consider_only]
    #stmapped.remove_species(remove_after_mapping)
    xsf_string = aenet.to_XSF(st_tmp, write_force_zero=True)
    xsf_string = "# total energy = {} eV\n\n".format(energies[i]) + xsf_string
    with open(os.path.join("aenet_xsf",
                           "perf_structure_stripped.{}.{}.xsf".format(myreplica, i)), "w") as fi:
        fi.write(xsf_string)


#perf_st = Structure.from_file("POSCAR")
#st = Structure.from_file("0/structure.0.vasp")

#stmapped = map2perflat(perf_st, st, vac_spaceholder)
#stmapped.to("POSCAR", "testmap.vasp")
