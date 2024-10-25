# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#               2024  Academeia Co., Ltd.
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

from pymatgen.core import Structure


"""
Adapted from pymatgen.io.xcrysden distributed under the MIT License
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""

def structure_to_XSF(structure: Structure, write_force_zero=False):
    """
    Returns a string with the structure in XSF format
    See http://www.xcrysden.org/doc/XSF.html
    """
    lines = []
    app = lines.append

    app("CRYSTAL")
    app("# Primitive lattice vectors in Angstrom")
    app("PRIMVEC")
    cell = structure.lattice.matrix
    for i in range(3):
        app(" %.14f %.14f %.14f" % tuple(cell[i]))

    cart_coords = structure.cart_coords
    app("# Cartesian coordinates in Angstrom.")
    app("PRIMCOORD")
    app(" %d 1" % len(cart_coords))
    species = structure.species
    site_properties = structure.site_properties
    if "forces" not in site_properties.keys():
        write_force_zero = True
    else:
        forces = site_properties["forces"]

    if write_force_zero:
        for a in range(len(cart_coords)):
            app(
                str(species[a])
                + " %20.14f %20.14f %20.14f" % tuple(cart_coords[a])
                + " 0.0 0.0 0.0"
            )
    else:
        for a in range(len(cart_coords)):
            app(
                str(species[a])
                + " %20.14f %20.14f %20.14f" % tuple(cart_coords[a])
                + " %20.14f %20.14f %20.14f" % tuple(forces[a])
            )

    return "\n".join(lines)


def structure_from_XSF(input_string: str):
    """
    Initialize a `Structure` object from a string with data in XSF format.

    Args:
        input_string: String with the structure in XSF format.
            See http://www.xcrysden.org/doc/XSF.html
        cls_: Structure class to be created. default: pymatgen structure

    """
    # CRYSTAL                                        see (1)
    # these are primitive lattice vectors (in Angstroms)
    # PRIMVEC
    #    0.0000000    2.7100000    2.7100000         see (2)
    #    2.7100000    0.0000000    2.7100000
    #    2.7100000    2.7100000    0.0000000

    # these are conventional lattice vectors (in Angstroms)
    # CONVVEC
    #    5.4200000    0.0000000    0.0000000         see (3)
    #    0.0000000    5.4200000    0.0000000
    #    0.0000000    0.0000000    5.4200000

    # these are atomic coordinates in a primitive unit cell  (in Angstroms)
    # PRIMCOORD
    # 2 1                                            see (4)
    # 16      0.0000000     0.0000000     0.0000000  see (5)
    # 30      1.3550000    -1.3550000    -1.3550000

    lattice, coords, species = [], [], []
    lines = input_string.splitlines()

    for i in range(len(lines)):
        if "PRIMVEC" in lines[i]:
            for j in range(i + 1, i + 4):
                lattice.append([float(c) for c in lines[j].split()])

        if "PRIMCOORD" in lines[i]:
            num_sites = int(lines[i + 1].split()[0])

            for j in range(i + 2, i + 2 + num_sites):
                tokens = lines[j].split()
                species.append(tokens[0])
                coords.append([float(j) for j in tokens[1:4]])
            break
    else:
        raise ValueError("Invalid XSF data")

    s = Structure(lattice, species, coords, coords_are_cartesian=True)
    return s
