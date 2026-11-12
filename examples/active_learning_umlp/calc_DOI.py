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

from pathlib import Path

from pymatgen.core import Structure
from abics.applications.latgas_abinitio_interface import naive_matcher

mapper = naive_matcher.naive_mapping
spinel_struct = Structure.from_file(Path(__file__).with_name("MgAl2O4.vasp"))

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
