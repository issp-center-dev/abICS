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

from mpi4py import MPI
from pymatgen import Structure
#from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
import os, sys
import numpy as np
import scipy.constants as constants
from abics.mc_mpi import RXParams
import naive_matcher

mapper = naive_matcher.naive_mapping

throwout = 0
tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
rxparams = RXParams.from_toml(tomlfile)
nreplicas = rxparams.nreplicas
comm = MPI.COMM_WORLD

spinel_struct = Structure.from_file("MgAl2O4.vasp")
spinel_struct.make_supercell([2,2,2])
Asite_struct = spinel_struct.copy()
Asite_struct.remove_species(["Al", "O"])


#matcher = StructureMatcher(
#    ltol=0.1,
#    primitive_cell=False,
#    allow_subset=True,
#    comparator=FrameworkComparator(),
#    ignored_species=["O"],
#)

asite_ids = spinel_struct.indices_from_symbol("Mg")
def calc_DOI(structure):
    mapping = mapper(spinel_struct, structure)
    x = 0
    species = structure.species
    species = [str(sp) for sp in species]
    for i in asite_ids:
        if species[mapping[i]] == "Al":
            x += 1
    x /= float(len(asite_ids))
    return x


myT = comm.Get_rank()
kTs = np.load(open("./kTs.npy", "rb"))

if myT < nreplicas:
    Trank_hist = np.load(open(str(myT) + "/Trank_hist.npy", "rb"))
    nstep = len(Trank_hist)
    os.chdir("Tseparate/{}".format(myT))
    DOI = []
    for i in range(throwout,nstep):
        DOI.append(calc_DOI(Structure.from_file("structure.{}.vasp".format(i))))
    DOI = np.array(DOI)
    np.savetxt("DOI.dat", DOI)
    DOImean = DOI.mean()
    os.chdir("..")
    myT_kelvin = kTs[myT] / constants.value(u"Boltzmann constant in eV/K")

for i in range(len(kTs)):
    comm.Barrier()
    if myT != i:
        continue
    outfi = open("DOI_T.dat", "a")
    outfi.write("{:6.1f} \t".format(myT_kelvin))
    outfi.write("{}\n".format(DOImean))
    outfi.close()

