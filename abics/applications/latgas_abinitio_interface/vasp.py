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
To deal with VASP
"""


from .base_solver import SolverBase
from collections import namedtuple
from pymatgen.io.vasp.inputs import Poscar, VaspInput
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
import numpy as np
import os.path


class VASPSolver(SolverBase):
    """
    VASP solver

    Attributes
    ----------
    path_to_solver : str
        Path to the solver
    input : VASPSolver.Input
        Input manager
    output : VASPSolver.Output
        Output manager
    """

    def __init__(self, path_to_solver):
        """
        Initialize the solver.

        Parameters
        ----------
        path_to_solver : str
            Path to the solver.
        """
        super(VASPSolver, self).__init__(path_to_solver)
        self.path_to_solver = path_to_solver
        self.input = VASPSolver.Input()
        self.output = VASPSolver.Output()

    def name(self):
        return "VASP"

    class Input(object):
        """
        Input manager for VASP

        Attributes
        ----------
        base_info : dict
            Common parameters
        pos_info : pymatgen.Structure
            Atom positions
        """

        def __init__(self):
            self.base_info = None
            self.pos_info = None

        def from_directory(self, base_input_dir):
            """

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            Returns
            -------
            base_vasp_input : VaspInput (defined in pymatgen)
                vasp input object
            """
            self.base_vasp_input = VaspInput.from_directory(base_input_dir)
            self.base_info = self.base_vasp_input.get("INCAR")
            return self.base_vasp_input

        def update_info_by_structure(self, structure):
            """
            Update information by structure file

            Parameters
            ----------
            structure : pymatgen.Structure
                Atomic structure
            """
            self.pos_info = Poscar(structure=structure, selective_dynamics=structure.site_properties["seldyn"])
            self.base_vasp_input.update({"POSCAR": self.pos_info})

        def update_info_from_files(self, workdir, rerun):
            """
            Update information from output files.

            Parameters
            ----------
            workdir: str
                Path to working directory.
            rerun: int
                Mode for rerunning (0: not rerun, 1: rerun (base file is replaced by new one))
            """
            if rerun == 1:
                info_dict = ["BASE", "POS"]
            elif rerun > 0:
                info_dict = ["POS"]

            # Add update procedure
            # 1. Check kind of the update information
            # 2. Replace the information
            for info in info_dict:
                if info == "BASE":
                    # Modify parameters
                    self.base_info = self.base_vasp_input.get("INCAR")
                    self.base_info.update({"IBRION": 3, "POTIM": 0.2})
                    self.base_vasp_input.update({"INCAR": self.base_info})
                elif info == "POS":
                    # Update positions
                    self.pos_info = Poscar.from_file(
                        os.path.join(workdir, "CONTCAR")
                    )
                    self.base_vasp_input.update({"POSCAR": self.pos_info})

        def write_input(self, output_dir):
            """
            Generate input files of the solver program.

            Parameters
            ----------
            workdir : str
                Path to working directory.
            """
            if self.base_info is None:
                raise AttributeError("Fail to set base_info.")
            self.base_vasp_input.write_input(output_dir=output_dir)

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
            self.drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
            self.queen = BorgQueen(self.drone)

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
            self.queen.serial_assimilate(workdir)
            results = self.queen.get_data()[-1]
            return Phys(np.float64(results.energy), results.structure)

    def solver_run_schemes(self):
        return ("mpi_spawn_ready",)
