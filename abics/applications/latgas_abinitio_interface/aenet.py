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

from __future__ import annotations

import os
import sys, shutil, io
from collections import namedtuple
import numpy as np
from pymatgen.core import Structure

from .base_solver import SolverBase
from .params import ALParams, DFTParams
from .util import structure_to_XSF, structure_from_XSF


class AenetSolver(SolverBase):
    """
    This class defines the aenet solver.
    """

    def __init__(
        self, path_to_solver: os.PathLike, ignore_species=None, run_scheme="subprocess"
    ):
        """
        Initialize the solver.

        Parameters
        ----------
        path_to_solver : str
                      Path to the solver.
        """
        super(AenetSolver, self).__init__(path_to_solver)
        self.path_to_solver = path_to_solver
        self.input = AenetSolver.Input(ignore_species, run_scheme)
        self.output = AenetSolver.Output()

    def name(self):
        return "aenet"

    class Input(object):
        def __init__(self, ignore_species: str | None, run_scheme="subprocess"):
            self.base_info = None
            self.pos_info = None
            self.ignore_species = ignore_species
            self.run_scheme = run_scheme

        def from_directory(self, base_input_dir: os.PathLike):
            """
            Initialize information from files in base_input_dir.

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            """

            # set information of base_input and pos_info from files in base_input_dir
            self.base_info = os.path.abspath(base_input_dir)
            # self.pos_info = open('{}/structure.xsf'.format(base_input_dir), 'r').read()

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
            self.pos_info = structure_to_XSF(structure)

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
            with open(os.path.join(output_dir, "structure.xsf"), "w") as f:
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
                    os.path.join(output_dir, "predict.in"),
                    os.path.join(output_dir, "structure.xsf"),
                    output_dir,
                ]
            elif self.run_scheme == "subprocess":
                return [
                    os.path.join(output_dir, "predict.in"),
                    os.path.join(output_dir, "structure.xsf"),
                ]
            else:
                raise RuntimeError("Invalid run_scheme: {}".format(self.run_scheme))

    class Output(object):
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
            with open(os.path.join(output_dir, "structure.xsf")) as f:
                structure = structure_from_XSF(f.read())
            with open(os.path.join(output_dir, "stdout")) as f:
                lines = f.read()
                fi_io = io.StringIO(lines)
                line = fi_io.readline()
            if "optimized" in lines:
                while "optimized" not in line:
                    line = fi_io.readline()
                for i in range(4):
                    fi_io.readline()
                for i in range(len(structure)):
                    xyz = [float(x) for x in fi_io.readline().split()[1:4]]
                    structure.replace(
                        i, structure[i].species, coords=xyz, coords_are_cartesian=True
                    )
            while "Total energy" not in line:
                line = fi_io.readline()
            energy = line.split()[-2]
            return Phys(np.float64(energy), structure)

    def solver_run_schemes(self):
        return ("subprocess", "mpi_spawn_ready")

    @classmethod
    def create(cls, params: ALParams | DFTParams):
        path = params.path
        ignore_species = params.ignore_species
        run_scheme = params.solver_run_scheme
        return cls(path, ignore_species, run_scheme)
