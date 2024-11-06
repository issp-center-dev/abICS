# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# abICS wrapper of NequIP solver
# Munehiro Kobayashi, Yusuke Konishi (Academeia Co., Ltd.) 2024
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
energy calculator using nequip python interface
"""

from __future__ import annotations

import os.path
from collections import namedtuple
import numpy as np
from pymatgen.core import Structure
import torch
from ase import Atoms
from nequip.data import AtomicDataDict, AtomicData
from nequip.utils import Config

from .base_solver import SolverBase, register_solver
from .params import ALParams, DFTParams


class NequipSolver(SolverBase):
    """
    Nequip solver

    Attributes
    ----------
    path_to_solver : str
        Path to the solver
    input : NequipSolver.Input
        Input manager
    output : NequipSolver.Output
        Output manager
    """

    def __init__(self, ignore_species):
        """
        Initialize the solver.

        """

        super(NequipSolver, self).__init__("")
        self.path_to_solver = self.calculate_energy
        self.input = NequipSolver.Input(ignore_species)
        self.output = NequipSolver.Output()

    def name(self):
        return "nequip"

    def calculate_energy(self, fi, output_dir):
        st = self.input.st
        symbols = [site.specie.symbol for site in st]
        positions = [site.coords for site in st]
        pbc = (True, True, True)
        cell = st.lattice.matrix
        atoms = Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell)

        atom_types = torch.tensor([self.input.element_list.index(atom) for atom in symbols], dtype=torch.long)

        data = AtomicData.from_ase(atoms, r_max=self.input.r_max)
        data[AtomicDataDict.ATOM_TYPE_KEY] = atom_types

        self.input.model.eval()
        with torch.no_grad():
            # Convert AtomicData to dictionary
            data_dict = data.to_dict()
            predicted = self.input.model(data_dict)

        # Get predicted energy
        ene = predicted['total_energy'].item()

        self.output.st = st
        self.output.ene = ene

    class Input(object):
        """
        Input manager for Mock

        Attributes
        ----------
        st : pymatgen.Structure
            structure
        """

        st: Structure

        def __init__(self, ignore_species=None):
            self.ignore_species = ignore_species
            # self.st = Structure()

        def from_directory(self, base_input_dir):
            """

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            """
            self.base_input_dir = base_input_dir
            self.model = torch.jit.load(os.path.join(base_input_dir, "deployed.pth"))
            yaml_file = os.path.join(base_input_dir, "input.yaml")
            yaml_dic = Config.from_file(yaml_file)
            self.element_list = yaml_dic["chemical_symbols"]
            self.r_max = yaml_dic["r_max"]

        def update_info_by_structure(self, structure):
            """
            Update information by structure file

            Parameters
            ----------
            structure : pymatgen.Structure
                Atomic structure
            """
            self.st = structure.copy()
            if self.ignore_species is not None:
                self.st.remove_species(self.ignore_species)

        def update_info_from_files(self, workdir, rerun):
            """
            Do nothing
            """
            pass

        def write_input(self, output_dir):
            """
            Generate input files of the solver program.

            Parameters
            ----------
            workdir : str
                Path to working directory.
            """
            if not os.path.exists(output_dir):
                import shutil

                shutil.copytree(self.base_input_dir, output_dir)

            # self.st.to("POSCAR", os.path.join(output_dir, "structure.vasp"))

        def cl_args(self, nprocs, nthreads, workdir):
            """
            Generate command line argument of the solver program.

            Parameters
            ----------
            nprocs : int
                The number of processes.
            nthreads : int
                The number of threads.
            workdir : str
                Path to the working directory.

            Returns
            -------
            args : list[str]
                Arguments of command
            """
            return [workdir]

    class Output(object):
        """
        Output manager.
        """

        def __init__(self):
            pass

        def get_results(self, workdir):
            """
            Get energy and structure obtained by the solver program.

            Parameters
            ----------
            workdir : str
                Path to the working directory.

            Returns
            -------
            phys : named_tuple("energy", "structure")
                Total energy and atomic structure.
                The energy is measured in the units of eV
                and coodinates is measured in the units of Angstrom.
            """
            Phys = namedtuple("PhysVaules", ("energy", "structure"))
            return Phys(self.ene, self.st)

    def solver_run_schemes(self):
        return ("function",)

    @classmethod
    def create(cls, params: ALParams | DFTParams):
        ignore_species = params.ignore_species
        return cls(ignore_species)
