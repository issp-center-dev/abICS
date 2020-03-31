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

from collections import namedtuple


class SolverBase(object):
    """
    Base class defining the common interface of solvers.

    Attributes
    ----------
    path_to_solver : str
        Path to solver program.
    input : SolverBase.Input
        Input manager.
    output : SolverBase.Output
        Output manager.
    """

    class Input(object):
        """
        Input manager.

        Attributes
        ----------
        base_info : Any
            Common parameter.
        pos_info : Any
            Information of position.
        """

        def __init__(self):
            self.base_info = None
            self.pos_info = None

        def update_info_by_structure(self, structure):
            """
            Update information by structure.

            Parameters
            ----------
            structure : pymatgen.Structure
                Atomic structure.

            """
            return None

        def update_info_from_files(self, workdir, rerun):
            """
            Update information by result files.

            Parameters
            ----------
            workdir : str
                Path to working directory.
            rerun : int
            """

            return None

        def write_input(self, workdir):
            """
            Generate input files of the solver program.

            Parameters
            ----------
            workdir : str
                Path to working directory.
            """
            return None

        def from_directory(self, base_input_dir):
            """
            Set information, base_input and pos_info, from files in the base_input_dir

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            """
            # set information of base_input and pos_info from files in base_input_dir
            self.base_info = {}
            self.pos_info = {}

        def cl_args(self, nprocs, nthreads, workdir):
            """
            Generate command line arguments of the solver program.

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
            return []

    class Output(object):
        """
        Output manager.
        """
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
            # Read results from files in workdir and calculate values
            phys = Phys(0.0, None)
            return phys

    def __init__(self, path_to_solver):
        """
        Initialize the solver.

        Parameters
        ----------
        solver_name : str
            Solver name.
        path_to_solver : str
            Path to the solver.
        """
        self.path_to_solver = path_to_solver
        self.input = SolverBase.Input
        self.output = SolverBase.Output

    def name(self):
        """
        Returns
        -------
        name : str
            Name of solver.
        """
        return "Base solver"

    def solver_run_schemes(self):
        """
        Returns
        -------
        schemes : tuple[str]
            Implemented runner schemes.
        """
        return ()
