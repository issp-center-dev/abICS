# from scipy.stats import norm
import numpy as np
import copy
from pymatgen import Lattice, Structure, Element, PeriodicSite
from pymatgen.io.vasp.inputs import Poscar
from mpi4py import MPI
from abics.applications.rhorhocy import rho_autocorr, rho_autocorr_fft


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


def rho_gauss(x, structure, species, sigma, comm=None):
    lattice = structure.lattice
    types_of_specie = [element.symbol for element in structure.types_of_specie]
    assert species in types_of_specie

    structure1 = structure.copy()

    not_species = copy.copy(types_of_specie)
    not_species.remove(species)

    structure1.remove_species(not_species)

    num_specie1 = structure1.num_sites
    numgrid = len(x)
    if comm:
        rho = np.zeros(numgrid)
        size = comm.Get_size()
        rank = comm.Get_rank()
        assert size < numgrid
        if numgrid % size != 0:
            raise ValueError("mesh size must be divisible by the communicator size")

        numgrid_loc = numgrid // size
        dist = lattice.get_all_distances(
            x[rank * numgrid_loc : (rank + 1) * numgrid_loc], structure1.frac_coords
        )
        rho[rank * numgrid_loc : (rank + 1) * numgrid_loc] = np.sum(
            np.exp(-np.power(dist / sigma, 2.0) / 2.0), axis=1
        )
        comm.Allreduce(MPI.IN_PLACE, rho, op=MPI.SUM)

    else:
        dist = lattice.get_all_distances(x, structure1.frac_coords)
        rho = np.sum(np.exp(-np.power(dist / sigma, 2.0) / 2.0), axis=1)

    rho /= np.sqrt(2.0 * np.pi) * sigma
    return rho


def write_rho_vasp(structure, species, sigma, mesh, filename):
    comm = MPI.COMM_WORLD
    size = comm.Get_rank()

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


if __name__ == "__main__":
    commworld = MPI.COMM_WORLD
    myrank = commworld.Get_rank()
    structure = Structure.from_file("POSCAR")
    species = "Y"
    sigma = 0.3
    mesh = [200, 200, 200]
    nsim = 24 * 24
    dr = 0.1
    rmax = 5.0
    if myrank == 0:
        write_rho_vasp(structure, species, sigma, mesh, "rhoY.vasp")
    structure.sort(key=lambda site: site.species_string)
    frac_coor_grid = np.array(
        np.meshgrid(range(mesh[0]), range(mesh[1]), range(mesh[2]), indexing="ij")
    ).T.reshape(-1, 3) / np.array(mesh)
    rho = rho_gauss(frac_coor_grid, structure, species, sigma)
    # autocorr = rho_autocorr(frac_coor_grid,rho,structure,dr,rmax,commworld,nsim)
    autocorr = rho_autocorr_fft(mesh, frac_coor_grid, rho, structure, dr, rmax)
    if myrank == 0:
        outfi = open("rhorho.dat", "w")
        outfi.write("\n".join([str(x) for x in autocorr]))
        outfi.close()
