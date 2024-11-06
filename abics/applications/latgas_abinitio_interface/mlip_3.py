# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# abICS wrapper of MLIP-3 solver
# Masashi Noda, Yusuke Konishi (Academeia Co., Ltd.) 2024
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

"""
Adapted from pymatgen.io.xcrysden distributed under the MIT License
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""

from __future__ import annotations

import io
import os
import shutil
import sys
from collections import namedtuple

import numpy as np
from pymatgen.core import Structure

from .base_solver import SolverBase
from .params import ALParams, DFTParams

def map_species_to_sequential_numbers(original_list):
    """
    Maps a list of species to sequential numbers, starting from 1.

    Parameters
    ----------
    original_list : list
        List of species.

    Returns
    -------
    list
        List of sequential numbers.
    """
    # Use a dictionary to map each unique element to a new number
    mapping = {}
    current_number = 1

    for item in original_list:
        if item not in mapping:
            mapping[item] = current_number
            current_number += 1

    # Map each element of the original list to the new number
    return [mapping[item] for item in original_list]

def to_CFG(structure: Structure, energy, write_force_zero=False):
    """
    Returns a string with the structure in CFG format
    CFG format is a format used in input of MLIP-3

    Parameters
    ----------
    structure : pymatgen.Structure
        Atomic structure
    energy : float
        Total energy
    write_force_zero : bool
        If True, the forces are written as zeros.
        If False, the forces are written as the forces in the structure object.
    
    Returns
    -------
    str
        String with the structure in CFG format
    """

    lines = []
    app = lines.append

    app("BEGIN_CFG")
    app(" Size")    

    cart_coords = structure.cart_coords
    app("%6d" % len(cart_coords))

    cell = structure.lattice.matrix
    app(" Supercell")
    for i in range(3):
        app("%16.6f%16.6f%16.6f" % tuple(cell[i]))

    species = structure.species
    mapped_species = map_species_to_sequential_numbers(species)

    site_properties = structure.site_properties
    if "forces" not in site_properties.keys():
        write_force_zero = True
    else:
        forces = site_properties["forces"]

    app(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz")
    if write_force_zero:
        for a in range(len(cart_coords)):
            app("%14d" % int(str(a+1))
                + "%5d" % mapped_species[a]
                + "%15.6f%14.6f%14.6f" % tuple(cart_coords[a])
                + "%13.6f%12.6f%12.6f" % tuple([0.0, 0.0, 0.0]) 
            )
    else:
        for a in range(len(cart_coords)):
            app("%14d" % int(str(a+1))
                + "%5d" % mapped_species[a]
                + "%15.6f%14.6f%14.6f" % tuple(cart_coords[a])
                + "%13.6f%12.6f%12.6f" % tuple(forces[a]) 
            )

    app(" Energy")
    app("%26.12f" % energy)
    app(" PlusStress:  xx          yy          zz          yz          xz          xy")
    app("%16.5f%12.5f%12.5f%12.5f%12.5f%12.5f" % tuple([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    app(" Feature   EFS_by       VASP")
    app("END_CFG")
    app("")
    app("")

    return "\n".join(lines)

def read_CFG(input_string: str):
    """
    Reads a string with the structure in CFG format and returns a dictionary

    Parameters
    ----------
    input_string : str
        String with the structure in CFG format

    Returns
    -------
    dict
        Dictionary with the structure in CFG format
    """
    cfg_dic = {}
    size = 0
    lines = input_string.split('\n')
    for i, line in enumerate(lines):
        if 'Size' in line:
            size = int(lines[i+1])
        if 'Supercell' in line:
            supercell = []
            for j in range(3):
                supercell.append([float(x) for x in lines[i+j+1].split()])
            cfg_dic['supercell'] = supercell
        if 'AtomData' in line:
            atom_data = []
            atom_type = []
            for j in range(size):
                atom_data.append([float(x) for x in lines[i+j+1].split()[2:5]])
                atom_type.append(int(lines[i+j+1].split()[1]))
            cfg_dic['atom_data'] = atom_data
            cfg_dic['atom_type'] = atom_type

    return cfg_dic

def from_CFG(input_string: str, species):
    """
    Returns a Structure object from a string with the structure in CFG format

    Parameters
    ----------
    input_string : str
        String with the structure in CFG format
    species : list
        List of species

    Returns
    -------
    pymatgen.Structure
        Atomic structure
    """
    cfg_dic = read_CFG(input_string)
    list_species = [species[i-1] for i in cfg_dic['atom_type']]
    s = Structure(cfg_dic['supercell'], list_species, cfg_dic['atom_data'], coords_are_cartesian=True)
    return s

# Need to mod for MLIP-3
class MLIP3Solver(SolverBase):
    """
    This class defines the MLIP-3 solver.
    """

    def __init__(
        self, path_to_solver: os.PathLike, ignore_species=None,
        run_scheme="subprocess"
    ):
        """
        Initialize the solver.

        Parameters
        ----------
        path_to_solver : str
                      Path to the solver.
        """
        super().__init__(path_to_solver)
        self.path_to_solver = path_to_solver
        self.species = None
        self.input = MLIP3Solver.Input(self, ignore_species, run_scheme)
        self.output = MLIP3Solver.Output(self)

    def name(self):
        return "mlip_3"

    class Input(object):
        def __init__(
            self, mlip3_solver, ignore_species: str | None, run_scheme="subprocess"
        ):
            self.mlip3_solver = mlip3_solver
            self.base_info = None
            self.pos_info = None
            self.pot_info = None
            self.ignore_species = ignore_species
            self.species = None
            self.run_scheme = run_scheme

        def from_directory(self, base_input_dir: os.PathLike):
            """
            Initialize information from files in base_input_dir.

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            """

            # set information of base_input and
            # pos_info from files in base_input_dir
            self.base_info = os.path.abspath(base_input_dir)
            # self.pos_info = open(
            #     '{}/structure.xsf'.format(base_input_dir), 'r'
            # ).read()

        def update_info_by_structure(self, structure: Structure):
            """
            Update information by atomic structure.

            Parameters
            ----------
            structure : pymatgen.Structure
                Atomic structure
            """
            if self.ignore_species is not None:
                structure = structure.copy()
                structure.remove_species(self.ignore_species)
            self.mlip3_solver.species = []
            seen = set()
            for specie in structure.species:
                if specie not in seen:
                    self.mlip3_solver.species.append(str(specie))
                    seen.add(specie)
            self.pos_info = to_CFG(structure, 0.0)

        def update_info_from_files(self, output_dir, rerun):
            """
            Do nothing.
            """
            print("rerun not implemented. Something has gone wrong")
            sys.exit(1)

        def write_input(self, output_dir: os.PathLike):
            """
            Generate input files of the solver program.

            Parameters
            ----------
            output_dir : os.PathLike
                Path to working directory.
            """
            # Write input files
            if self.base_info is None:
                raise AttributeError("Fail to set base_info.")
            os.makedirs(output_dir, exist_ok=True)
            for fname in os.listdir(self.base_info):
                shutil.copy(os.path.join(self.base_info, fname), output_dir)
            with open(os.path.join(output_dir, "structure.cfg"), "w") as f:
                f.write(self.pos_info)

        def cl_args(self, nprocs, nthreads, output_dir):
            """
            Generate command line arguments of the solver program.

            Parameters
            ----------
            nprocs : int
                The number of processes.
            nthreads : int
                The number of threads.
            output_dir : str
                Path to the working directory.

            Returns
            -------
            args : list[str]
                Arguments of command
            """
            # Specify command line arguments
            if self.run_scheme == "mpi_spawn_ready":
                return [
                    "calculate_efs",
                    os.path.join(output_dir, "pot.almtp"),
                    os.path.join(output_dir, "structure.cfg"),
                    output_dir,
                ]
            elif self.run_scheme == "subprocess":
                return [
                    "calculate_efs",
                    os.path.join(output_dir, "pot.almtp"),
                    os.path.join(output_dir, "structure.cfg"),
                ]

    class Output(object):
        def __init__(self, mlip3_solver):
            self.mlip3_solver = mlip3_solver

        def get_results(self, output_dir):
            """
            Get energy and structure obtained by the solver program.

            Parameters
            ----------
            output_dir: str
                Path to the working directory.

            Returns
            -------
            phys : named_tuple("energy", "structure")
                Total energy and atomic structure.
                The energy is measured in the units of eV
                and coordinates is measured in the units of Angstrom.

            """
            # Read results from files in output_dir and calculate values
            Phys = namedtuple("PhysValues", ("energy", "structure"))
            with open(os.path.join(output_dir, "structure.cfg")) as f:
                lines = f.read()
            structure = from_CFG(lines, self.mlip3_solver.species)
            fi_io = io.StringIO(lines)
            line = fi_io.readline()
            #if "optimized" in lines:
            #    while "optimized" not in line:
            #        line = fi_io.readline()
            #    for i in range(4):
            #        fi_io.readline()
            #    for i in range(len(structure)):
            #        xyz = [float(x) for x in fi_io.readline().split()[1:4]]
            #        structure.replace(
            #            i, structure[i].species, coords=xyz,
            #            coords_are_cartesian=True
            #        )
            while "Energy" not in line:
                line = fi_io.readline()
            line = fi_io.readline()
            energy = line
            return Phys(np.float64(energy), structure)

    def solver_run_schemes(self):
        return ("subprocess", "mpi_spawn_ready")

    @classmethod
    def create(cls, params: ALParams | DFTParams):
        path = params.path
        ignore_species = params.ignore_species
        run_scheme = params.solver_run_scheme
        return cls(path, ignore_species, run_scheme)
