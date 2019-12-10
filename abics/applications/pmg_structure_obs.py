import numpy as np
import copy
from pymatgen import Lattice, Structure, Element, PeriodicSite
from pymatgen.io.vasp.inputs import Poscar


def g_r(structure, specie1, specie2, grid_1D):
    X = grid_1D.x
    dr = grid_1D.dx
    assert grid_1D.x[0] == dr
    lattice = structure.lattice
    types_of_specie = [element.symbol for element in structure.types_of_specie]
    assert specie1 in types_of_specie
    assert specie2 in types_of_specie

    structure1 = structure.copy()
    structure2 = structure.copy()

    # print(specie1)
    not_specie1 = copy.copy(types_of_specie)
    not_specie1.remove(specie1)
    not_specie2 = types_of_specie
    not_specie2.remove(specie2)
    # print(not_specie1,not_specie2)

    structure1.remove_species(not_specie1)
    structure2.remove_species(not_specie2)

    num_specie1 = structure1.num_sites
    num_specie2 = structure2.num_sites

    dist = lattice.get_all_distances(structure1.frac_coords, structure2.frac_coords)
    dist_bin = np.around(dist / dr)
    # print(num_specie1,num_specie2)
    g = np.zeros(len(X))
    for i in range(len(X)):
        g[i] = np.count_nonzero(dist_bin == i + 1)

    if specie1 == specie2:
        pref = lattice.volume / (
            4.0 * np.pi * X * X * dr * num_specie1 * (num_specie1 - 1)
        )
    else:
        pref = lattice.volume / (4.0 * np.pi * X * X * dr * num_specie1 * num_specie2)

    g *= pref
    return g


def isotropic_gauss_nd(x, xcenter, sigma):
    return np.prod(
        1.0
        / (np.sqrt(2.0 * np.pi * sigma))
        * np.exp(-np.power((x - xcenter) / sigma, 2.0) / 2.0)
    )


def gauss_1d(x, center, sigma):
    return (
        1.0
        / (np.sqrt(2.0 * np.pi) * sigma)
        * np.exp(-np.power((x - xcenter) / sigma, 2.0) / 2.0)
    )


def write_rho_vasp(structure, species, sigma, mesh, filename):
    frac_coor_grid = np.array(
        np.meshgrid(range(mesh[0]), range(mesh[1]), range(mesh[2]), indexing="ij")
    ).T.reshape(-1, 3) / np.array(mesh)

    # np.array(list(np.ndindex(mesh[0],mesh[1],mesh[2])))/np.array(mesh)
    structure1 = structure.copy()
    structure1.sort(key=lambda site: site.species_string)
    poscar = Poscar(structure1)
    outfi = open(filename, "w")
    print(frac_coor_grid)
    rho = rho_gauss(frac_coor_grid, structure, species, sigma)
    outfi.write("rho_vasp\n")
    outfi.write("\t1.00\n")
    latticevecs = structure.lattice.matrix
    outfi.write("\t".join([str(l) for l in latticevecs[0]]) + "\n")
    outfi.write("\t".join([str(l) for l in latticevecs[1]]) + "\n")
    outfi.write("\t".join([str(l) for l in latticevecs[2]]) + "\n")
    outfi.write("\t".join(poscar.site_symbols) + "\n")
    outfi.write("\t".join([str(l) for l in poscar.natoms]) + "\n")
    outfi.write("Direct\n")
    for site in poscar.structure:
        coords = site.frac_coords
        outfi.write("\t".join([str(l) for l in coords]) + "\n")
    outfi.write("\n")
    outfi.write("\t".join(str(l) for l in mesh) + "\n")
    for i in range(len(rho) // 5 + 1):
        outfi.write(" ".join(str(x) for x in rho[i * 5 : (i + 1) * 5]) + "\n")
