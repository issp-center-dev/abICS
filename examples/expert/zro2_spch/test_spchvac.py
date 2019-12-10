import numpy as np
import random as rand
import sys, os
import copy
import pickle

# from mpi4py import MPI

from pymatgen import Lattice, Structure, Element
from pymatgen.io.vasp import Poscar, VaspInput
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen


if __name__ == "__main__":

    # prepare spch_config
    ltol = 0.01
    N_ovac = 0.05
    base_structure = Structure.from_file(
        os.path.join(os.path.dirname(__file__), "POSCAR")
    )  # .get_primitive_structure(tolerance=0.001)
    base_Osite = base_structure.copy()
    base_Osite.remove_species(["Pt", "Zr"])

    Osites = base_structure.indices_from_symbol("O")
    N_Vsites = int(N_ovac * len(Osites))
    Vacsites = rand.sample(Osites, N_Vsites)
    Vacsites.sort()
    Vacsites.reverse()
    print(Vacsites)
    structure = base_structure.copy()
    for site in Vacsites:
        structure.pop(site)

    structure_Osite = structure.copy()
    structure_Osite.remove_species(["Pt", "Zr"])

    matcher_site = StructureMatcher(
        ltol=ltol,
        primitive_cell=False,
        allow_subset=True,
        comparator=FrameworkComparator(),
        ignored_species=["Pt"],
    )

    # Figure out where the vacancies are by comparing with base_structure
    filled_sites = matcher_site.get_mapping(base_Osite, structure_Osite)

    print(filled_sites)

    vac_struct = base_Osite.copy()
    vac_struct.remove_sites(filled_sites)
    print(vac_struct)
    # assert len(vac_struct) == 10

    # choose one vacancy and one O atom randomly and flip
    vac_flip = rand.choice(vac_struct)
    Osites = structure.indices_from_symbol("O")
    O_flip = rand.choice(Osites)

    del structure[O_flip]
    print(vac_flip.frac_coords)
    structure.append("O", vac_flip.frac_coords)
    print(structure)
