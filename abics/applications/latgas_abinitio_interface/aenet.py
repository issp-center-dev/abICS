from .base_solver import SolverBase
from collections import namedtuple
import numpy as np
from pymatgen import Structure
import os
import sys,shutil,io

class aenetSolver(SolverBase):
    """
    This class defines the aenet solver.
    """

    def __init__(self, path_to_solver):
        """
        Initialize the solver.

        Parameters
        ----------
        path_to_solver : str
                      Path to the solver.
        """
        super(aenetSolver, self).__init__(path_to_solver)
        self.path_to_solver = path_to_solver
        self.input = aenetSolver.Input()
        self.output = aenetSolver.Output()

    def name(self):
        return "aenet"

    class Input(object):
        def __init__(self):
            self.base_info = None
            self.pos_info = None

        def from_directory(self, base_input_dir):
            """

            Parameters
            ----------
            base_input_dir

            Returns
            -------

            """
            # set information of base_input and pos_info from files in base_input_dir
            self.base_info = os.path.abspath(base_input_dir)
            self.pos_info = open('{}/structure.xsf'.format(base_input_dir), 'r').read()

        def update_info_by_structure(self, structure):
            """

            Parameters
            ----------
            structure

            Returns
            -------

            """
            self.pos_info = structure.to('XSF')

        def update_info_from_files(self, output_dir, rerun):
            """

            Parameters
            ----------
            output_dir
            rerun

            Returns
            -------

            """
            print('rerun not implemented. Something has gone wrong')
            sys.exit(1)
            
        def write_input(self, output_dir):
            """

            Parameters
            ----------
            output_dir

            Returns
            -------

            """
            # Write input files
            if self.base_info is None:
                raise AttributeError("Fail to set base_info.")
            os.makedirs(output_dir, exist_ok=True)
            for fname in os.listdir(self.base_info):
                shutil.copy('{}/{}'.format(self.base_info, fname), output_dir)
            with open('{}/structure.xsf'.format(output_dir), 'w') as f:
                f.write(self.pos_info)

        def cl_args(self, nprocs, nthreads, output_dir):
            """

            Parameters
            ----------
            nprocs
            nthreads
            output_dir

            Returns
            -------

            """
            # Specify command line arguments
            return ['{}/predict.in'.format(output_dir),]

    class Output(object):

        def get_results(self, output_dir):
            """

            Parameters
            ----------
            output_dir

            Returns
            -------

            """
            # Read results from files in output_dir and calculate values
            Phys = namedtuple("PhysVaules", ("energy", "structure"))
            structure = Structure.from_file('{}/structure.xsf'.format(output_dir))
            with open('{}/stdout'.format(output_dir)) as f:
                lines = f.read()
                fi_io = io.StringIO(lines)
                line = fi_io.readline()
            if 'optimized' in lines:
                while 'optimized' not in line:
                    line = fi_io.readline()
                for i in range(4):
                    fi_io.readline()
                for i in range(len(structure)):
                    xyz = [float(x) for x in fi_io.readline().split()[1:4]]
                    structure.replace(i, structure[i].species, coords = xyz,
                                      coords_are_cartesian = True)
            while 'Total energy' not in line:
                line = fi_io.readline()
            energy = line.split()[-2]
            return Phys(np.float64(energy), structure)

    def solver_run_schemes(self):
        return ('subprocess')

