# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# abICS wrapper of CHGNet
# Yusuke Konishi (Academeia Co., Ltd.) 2025
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
energy calculator using sevennet python interface
"""

from __future__ import annotations

import os.path
from collections import namedtuple
from pymatgen.core import Structure
import torch
from ase import Atoms
from ase.calculators.calculator import Calculator
from nequip.utils import Config
from chgnet.model.model import CHGNet
from chgnet.model import StructOptimizer

from .base_solver import SolverBase, register_solver
from .params import ALParams, DFTParams

class CHGNetSolver(SolverBase):
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

        super(CHGNetSolver, self).__init__("")
        self.path_to_solver = self.calculate_energy
        self.input = CHGNetSolver.Input(ignore_species, use_pretrained, relax, fmax, device)
        self.output = CHGNetSolver.Output()

    def name(self):
        return "chgnet"

    def calculate_energy(self, fi, output_dir):
        st = self.input.st

        if self.input.relax:
            result = self.input.calculator.relax(atoms=st, fmax=self.input.fmax)
            # Get predicted energy
            ene = result["trajectory"].energies[-1]
        else:
            comp = st.composition.as_dict()
            num_atoms = sum([comp[key] for key in comp.keys()])
            ene = self.input.calculator.predict_structure(st)["e"]*num_atoms

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

        def __init__(self, ignore_species=None, use_pretrained=True, relax=True, fmax=0.1, device="cpu"):
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
            model_file = os.path.join(base_input_dir, "deployed.pth.tar")
            if not self.use_pretrained:
                if self.relax:
                    model = CHGNet.from_file(model_file)
                    self.calculator = StructOptimizer(model=model, optimizer_class="BFGS")
                else:
                    self.calculator = CHGNet.from_file(model_file)
            else:
                if self.relax:
                    model = CHGNet.load(use_device=self.device)
                    self.calculator = StructOptimizer(model=model, optimizer_class="BFGS")
                else:
                    self.calculator = CHGNet.load(use_device=self.device)

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
