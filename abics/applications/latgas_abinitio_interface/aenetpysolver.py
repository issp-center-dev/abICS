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
energy calculator using aenet python interface
"""

import numpy as np

from .base_solver import SolverBase
from collections import namedtuple
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
import os.path

try:
    import ase.io
except ImportError:
    sys.stderr.write("Error: unable to import ASE.")
    sys.exit()

from aenet.ase_calculator import ANNCalculator

import time

class AenetPySolver(SolverBase):
    """
    Aenetpy solver

    Attributes
    ----------
    path_to_solver : str
        Path to the solver
    input : AenetPySolver.Input
        Input manager
    output : AenetPySolver.Output
        Output manager
    """

    def __init__(self, potentials, ignore_species):
        """
        Initialize the solver.

        """
        super(AenetPySolver, self).__init__("")
        self.path_to_solver = self.calculate_energy
        self.input = AenetPySolver.Input(ignore_species)
        self.output = AenetPySolver.Output()
        print(potentials)
        self.calc = ANNCalculator(potentials)
        

    def name(self):
        return "aenetPy"

    def calculate_energy(self, fi, output_dir):
        times = []
        times.append(time.time())
        st = self.input.st
        times.append(time.time())
        atoms = AseAtomsAdaptor.get_atoms(st)
        times.append(time.time())
        atoms.set_calculator(self.calc)
        times.append(time.time())
        ene = atoms.get_potential_energy()
        times.append(time.time())
        self.output.st = st
        self.output.ene = ene
        times.append(time.time())
        dt = []
        for i in range(1,len(times)):
            dt.append(times[i] - times[i-1])
        print("dt",dt)

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
            os.makedirs(output_dir, exist_ok=True)
            #self.st.to("POSCAR", os.path.join(output_dir, "structure.vasp"))

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
