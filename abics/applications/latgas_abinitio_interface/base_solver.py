from collections import namedtuple


class SolverBase(object):
    """
    This class defines the common interface of solvers.
    """

    class Input(object):
        def __init__(self):
            self.base_info = None
            self.pos_info = None

        def update_info_by_structure(self, structure):

            return None

        def update_info_from_files(self, output_dir, rerun):
            if rerun == 1:
                info_dict = ["BASE", "POS"]
            elif rerun > 0:
                info_dict = ["BASE"]
            # 1. Check kind of the update information
            # 2. Replace the information by reading files in output_dir
            # for info in info_dict:
            #    if key == "BASE":
            #        #Modify parameters
            #        self.base_info = update_base_info
            #    elif key == "POS":
            #        #Modify parameters
            #       self.pos_info = update_pos_info
            return None

        def write_input(self, output_dir):
            # Write input files
            return None

        def from_directory(self, base_input_dir):
            # set information of base_input and pos_info from files in base_input_dir
            base_info = {}
            pos_info = {}
            self.base_info = {}
            self.pos_info = {}

        def cl_args(self, nprocs, nthreads, workdir):
            return []
            
    class Output(object):

        def get_results(self, output_dir):
            Phys = namedtuple("PhysVaules", ("energy", "structure"))
            # Read results from files in output_dir and calculate values
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
        return "Base solver"

    def solver_run_schemes(self):
        return ()
