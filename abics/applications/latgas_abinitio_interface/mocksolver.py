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
Mock for energy calculator
"""


from .base_solver import SolverBase
from collections import namedtuple
from pymatgen.core import Structure
import os.path


class MockSolver(SolverBase):
    """
    Mock solver

    Attributes
    ----------
    path_to_solver : str
        Path to the solver
    input : MockSolver.Input
        Input manager
    output : MockSolver.Output
        Output manager
    """

    def __init__(self):
        """
        Initialize the solver.

        """
        super(MockSolver, self).__init__("")
        self.path_to_solver = self.calc_energy
        self.input = MockSolver.Input()
        self.output = MockSolver.Output()

    def name(self):
        return "Mock"

    def calc_energy(self, fi, output_dir):
        st = Structure.from_file(os.path.join(output_dir, "structure.vasp"))
        ene = 0.0
        for sp in st.species:
            st_local = st.copy()
            st_local.remove_species(filter(lambda s: s != sp, st.species))
            dm = st_local.distance_matrix
            n = len(st_local)
            for i in range(n):
                for j in range(i+1,n):
                    ene += dm[i,j] ** 2
        with open(os.path.join(output_dir, "energy.dat"), "w") as f:
            f.write('{:.15f}\n'.format(ene))

    class Input(object):
        """
        Input manager for Mock

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
            self.st.to("POSCAR", os.path.join(output_dir, "structure.vasp"))

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
            Phys = namedtuple("PhysVaules", ("energy", "structure"))
            with open(os.path.join(workdir, "energy.dat")) as f:
                ene = float(f.read())
            st = Structure.from_file(os.path.join(workdir, "structure.vasp"))
            return Phys(ene, st)

    def solver_run_schemes(self):
        return ("function",)
