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

from collections import namedtuple

from pymatgen.core.structure import Structure

from .params import ALParams, DFTParams

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

        def __init__(self) -> None:
            self.base_info = None
            self.pos_info = None

        def update_info_by_structure(self, structure: Structure) -> None:
            """
            Update information by structure.

            Parameters
            ----------
            structure : pymatgen.Structure
                Atomic structure.

            """
            return None

        def update_info_from_files(self, workdir, rerun) -> None:
            """
            Update information by result files.

            Parameters
            ----------
            workdir : str
                Path to working directory.
            rerun : int
            """

            return None

        def write_input(self, workdir) -> None:
            """
            Generate input files of the solver program.

            Parameters
            ----------
            workdir : str
                Path to working directory.
            """
            return None

        def from_directory(self, base_input_dir) -> None:
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

        def cl_args(self, nprocs, nthreads, workdir) -> list[str]:
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
            ret: list[str] = []
            return ret

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
            Phys = namedtuple("PhysValues", ("energy", "structure"))
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

    @classmethod
    def create(cls, params: ALParams | DFTParams) -> SolverBase:
        """
        Create solver instance.

        Parameters
        ----------
        params : ALParams or DFTParams
            Parameters.

        Returns
        -------
        solver : SolverBase
            Solver instance.
        """
        raise NotImplementedError()


__solver_table = {}

def register_solver(solver_name: str, solver_class: str, solver_module: str) -> None:
    """
    Register solver class.

    Parameters
    ----------
    solver_name : str
        Solver name (case insensible).
    solver_class : str
        Solver class, which should be a subclass of SolverBase.
    solver_module : str
        Module name including the solver class.
    """

    __solver_table[solver_name.lower()] = (solver_class, solver_module)


def create_solver(solver_name, params: ALParams | DFTParams) -> SolverBase:
    """
    Create solver instance.

    Parameters
    ----------
    solver_name : str
        Solver name (case insensible).
    params : ALParams or DFTParams
        Parameters.

    Returns
    -------
    solver : SolverBase
        Solver instance.
    """
    sn = solver_name.lower()
    if sn not in __solver_table:
        raise ValueError(f"Unknown solver: {solver_name}")

    import importlib
    solver_class_name, solver_module = __solver_table[sn]
    mod = importlib.import_module(solver_module)
    solver_class = getattr(mod, solver_class_name)
    if SolverBase not in solver_class.mro():
        raise TypeError("solver_class must be a subclass of SolverBase")
    return solver_class.create(params)