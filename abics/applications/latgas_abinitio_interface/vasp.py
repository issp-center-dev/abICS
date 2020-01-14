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
            set information of base_input and pos_info from files in base_input_dir
            """
            self.base_vasp_input = VaspInput.from_directory(base_input_dir)
            self.base_info = self.base_vasp_input.get("INCAR")
            return self.base_vasp_input

        def update_info_by_structure(self, structure, seldyn_arr=None):
            self.pos_info = Poscar(structure=structure, selective_dynamics=seldyn_arr)
            self.base_vasp_input.update({"POSCAR": self.pos_info})

        def update_info_from_files(self, workdir, rerun):
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

        def write_input(self, workdir):
            """
            Generate input files of the solver program.

            Parameters
            ----------
            workdir : str
                Path to working directory.
            """
            if self.base_info is None:
                raise AttributeError("Fail to set base_info.")
            self.base_vasp_input.write_input(output_dir=workdir)

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
