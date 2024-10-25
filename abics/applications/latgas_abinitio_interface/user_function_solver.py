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

"""
User-defined function as energy calculator
"""

from __future__ import annotations

from collections import namedtuple
import os
import sys
import importlib

from pymatgen.core import Structure

from .base_solver import SolverBase
from .params import ALParams, DFTParams


class UserFunctionSolver(SolverBase):
    """
    UserFunction solver

    Attributes
    ----------
    path_to_solver : function(st) -> float
    input : SimpleFunctionSolver.Input
        Input manager
    output : SimpleFunctionSolver.Output
        Output manager
    """

    def __init__(self, function = None):
        """
        Initialize the solver.

        """
        super(UserFunctionSolver, self).__init__("")
        self.input = UserFunctionSolver.Input()
        self.output = UserFunctionSolver.Output()
        self.path_to_solver = self.__call__
        self.__fn = function

    def __call__(self, fi, output_dir):
        st = self.input.st
        if self.__fn is None:
            raise ValueError("Function is not set")
        ene = self.__fn(st)
        with open(os.path.join(output_dir, "energy.dat"), "w") as f:
            f.write("{:.15f}\n".format(ene))

    def name(self):
        return "UserFunction"

    def set_function(self, fn):
        self.__fn = fn

    class Input(object):
        """
        Input manager for SimpleFunction

        Attributes
        ----------
        st : pymatgen.Structure
            structure
        """

        st: Structure

        def __init__(self):
            # self.st = Structure()
            pass

        def from_directory(self, base_input_dir):
            """

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            """
            pass

        def update_info_by_structure(self, structure):
            """
            Update information by structure file

            Parameters
            ----------
            structure : pymatgen.Structure
                Atomic structure
            """
            self.st = structure

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
            os.makedirs(output_dir, exist_ok=True)
            self.st.to(
                fmt="POSCAR", filename=os.path.join(output_dir, "structure.vasp")
            )

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
            # Read results from files in output_dir and calculate values
            Phys = namedtuple("PhysValues", ("energy", "structure"))
            with open(os.path.join(workdir, "energy.dat")) as f:
                ene = float(f.read())
            st = Structure.from_file(os.path.join(workdir, "structure.vasp"))
            return Phys(ene, st)

    def solver_run_schemes(self):
        return ("function",)

    @classmethod
    def create(cls, params: ALParams | DFTParams):
        fn = None
        if params.function_module:
            sys.path.append(os.getcwd())
            pkg = params.function_module.split('.')
            modname = '.'.join(pkg[:-1])
            funcname = pkg[-1]
            try:
                mod = importlib.import_module(modname)
                fn = getattr(mod, funcname)
            except ImportError:
                raise ImportError(f"Cannot import module {modname}")
        return cls(fn)
