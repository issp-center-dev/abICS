# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# abICS wrapper of Mace
# Masashi Noda, Yusuke Konishi (Academeia Co., Ltd.) 2025
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
energy calculator using mace python interface
"""

from __future__ import annotations

import os.path
from collections import namedtuple
import numpy as np
from pymatgen.core import Structure
import torch
from ase import Atoms
from ase.optimize import BFGS
from nequip.data import AtomicDataDict, AtomicData
from nequip.utils import Config
from mace.calculators import mace_mp, MACECalculator

from .base_solver import SolverBase, register_solver
from .params import ALParams, DFTParams


class MaceSolver(SolverBase):
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

    def __init__(self, ignore_species, use_pretrained, relax, fmax, device):
        """
        Initialize the solver.

        """

        super(MaceSolver, self).__init__("")
        self.path_to_solver = self.calculate_energy
        self.input = MaceSolver.Input(ignore_species, use_pretrained, relax, fmax, device)
        self.output = MaceSolver.Output()

    def name(self):
        return "mace"

    def calculate_energy(self, fi, output_dir):
        st = self.input.st
        symbols = [site.specie.symbol for site in st]
        positions = [site.coords for site in st]
        pbc = (True, True, True)
        cell = st.lattice.matrix
        atoms = Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell)

        atoms.calc = self.input.calculator

        if self.input.relax:
            opt = BFGS(atoms)
            opt.run(fmax=self.input.fmax)

        # Get predicted energy
        ene = atoms.get_total_energy()

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

        def __init__(self, ignore_species=None, use_pretrained=False, relax=False, fmax=0.1, device="cpu"):
            self.ignore_species = ignore_species
            self.use_pretrained = use_pretrained
            self.relax = relax
            self.fmax = fmax
            self.device = device
            # self.st = Structure()

        def from_directory(self, base_input_dir):
            """

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            """
            self.base_input_dir = base_input_dir
            if not self.use_pretrained:
                self.model = os.path.join(base_input_dir, "deployed.model")
                self.calculator = MACECalculator(model_paths=self.model, device=self.device, default_dtype="float64")
            else:
                self.calculator = mace_mp(device=self.device, default_dtype="float64")

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
        use_pretrained = params.use_pretrained
        relax = params.relax
        fmax = params.fmax
        device = params.device
        return cls(ignore_species, use_pretrained, relax, fmax, device)
