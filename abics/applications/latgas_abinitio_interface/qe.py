from collections import namedtuple
import xml.etree.ElementTree as ET
import os.path

import numpy as np
import scipy.constants as spc
from pymatgen import Structure
from qe_tools.parsers import PwInputFile

from .base_solver import SolverBase
from ...util import expand_path


hartree2eV = spc.value("Hartree energy in eV")
Bohr2AA = spc.value('Bohr radius') * 1e10

def load_QE_result(path):
    return ComputedEntry(composition=Composition(comp), energy=ene)


class QESolver(SolverBase):
    """
    This class defines the QE solver.
    """

    def __init__(self, path_to_solver):
        """
        Initialize the solver.

        Parameters
        ----------
        path_to_solver : str
                      Path to the solver.
        """
        super(QESolver, self).__init__(path_to_solver)
        self.path_to_solver = path_to_solver
        self.input = QESolver.Input()
        self.output = QESolver.Output("pwscf")

    def name(self):
        return "QuantumESPRESSO"

    class Input(object):
        def __init__(self):
            self.pwi = None
            self.datadir = "pwscf.save"
            self.filetocheck = "data-file-schema.xml"

        def cleanup(self, rundir):
            checkfile = os.path.join(rundir, self.datadir, self.filetocheck)
            if os.path.exists(checkfile):
                os.remove(checkfile)

        def check_finished(self, rundir):
            f = os.path.join(rundir, self.datadir, self.filetocheck)
            if not os.path.exists(f):
                return False
            try:
                ET.parse(f)
            except ET.ParseError:
                return False
            return True

        def from_directory(self, base_input_dir):
            # set information of base_input and pos_info from files in base_input_dir

            # self.base_vasp_input = VaspInput.from_directory(base_input_dir)
            # self.base_info = self.base_vasp_input.get("INCAR")

            self.pwi = PwInputFile(os.path.join(os.getcwd(), base_input_dir, "scf.in"))
            self.pwi.namelists["CONTROL"]["prefix"] = "pwscf"
            self.pwi.namelists["CONTROL"]["pseudo_dir"] = expand_path(
                self.pwi.namelists["CONTROL"]["pseudo_dir"], os.getcwd()
            )

        def update_info_by_structure(self, structure, seldyn_arr=None):
            A = structure.lattice.matrix
            self.pwi.cell_parameters = {"units": "alat", "cell": A}

            nat = len(structure.sites)
            self.pwi.namelists["SYSTEM"]["ibrav"] = 0
            self.pwi.namelists["SYSTEM"]["a"] = 1.0
            self.pwi.namelists["SYSTEM"]["ntyp"] = len(structure.composition.elements)
            self.pwi.namelists["SYSTEM"]["nat"] = nat
            self.pwi.atomic_positions["names"] = [""] * nat
            self.pwi.atomic_positions["positions"] = [[0.0, 0.0, 0.0]] * nat
            self.pwi.atomic_positions["fixed_coords"] = [[False, False, False]] * nat

            for i, site in enumerate(structure.sites):
                self.pwi.atomic_positions["names"][i] = site.specie
                self.pwi.atomic_positions["positions"][i] = [site.a, site.b, site.c]
                self.pwi.atomic_positions["fixed_coords"][i] = [False, False, False]

        def update_info_from_files(self, output_dir, rerun):
            # TODO
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
                        os.path.join(output_dir, "CONTCAR")
                    )
                    self.base_vasp_input.update({"POSCAR": self.pos_info})
            """

        def write_input(self, output_dir):
            self.pwi.namelists["CONTROL"]["outdir"] = output_dir
            os.makedirs(output_dir, exist_ok=True)
            with open(os.path.join(output_dir, "scf.in"), "w") as f:
                for section in self.pwi.namelists.keys():
                    f.write("&{}\n".format(section))
                    for k, v in self.pwi.namelists[section].items():
                        if isinstance(v, str):
                            f.write('  {} = "{}"\n'.format(k, v))
                        elif isinstance(v, bool):
                            f.write(
                                "  {} = {}\n".format(k, ".true." if v else ".false.")
                            )
                        else:
                            f.write("  {} = {}\n".format(k, v))
                    f.write("/\n")

                f.write(
                    "CELL_PARAMETERS {}\n".format(self.pwi.cell_parameters["units"])
                )
                for row in self.pwi.cell_parameters["cell"]:
                    for elem in row:
                        f.write(" {}".format(elem))
                    f.write("\n")

                f.write("ATOMIC_SPECIES\n")
                for name, mass, pfile in zip(
                    self.pwi.atomic_species["names"],
                    self.pwi.atomic_species["masses"],
                    self.pwi.atomic_species["pseudo_file_names"],
                ):
                    f.write(" {} {} {}\n".format(name, mass, pfile))

                f.write(
                    "ATOMIC_POSITIONS {}\n".format(self.pwi.atomic_positions["units"])
                )
                for name, pos, if_pos in zip(
                    self.pwi.atomic_positions["names"],
                    self.pwi.atomic_positions["positions"],
                    self.pwi.atomic_positions["fixed_coords"],
                ):
                    f.write(" {}".format(name))
                    for p in pos:
                        f.write(" {}".format(p))
                    for p in if_pos:
                        if p:
                            f.write(" 0")
                        else:
                            f.write(" 1")
                    f.write("\n")

                f.write("K_POINTS {}\n".format(self.pwi.k_points["type"]))
                for kp in self.pwi.k_points["points"]:
                    f.write(" {}".format(kp))
                for offset in self.pwi.k_points["offset"]:
                    f.write(" {}".format(int(2 * offset)))
                f.write("\n")

        def file_for_check_finished(self, output_dir):
            return os.path.join(output_dir, self.datadir, self.filetocheck)

        def cl_args(self, nprocs, nthreads, output_dir):
            # Specify command line arguments
            return ["-in", os.path.join(output_dir, "scf.in")]

    class Output(object):
        def __init__(self, prefix):
            pass

        def get_results(self, output_dir):
            # Read results from files in output_dir and calculate values
            tree = ET.parse(
                os.path.join(output_dir, "pwscf.save", "data-file-schema.xml")
            )
            root = tree.getroot()
            A = np.zeros((3, 3))
            child = root.find("input").find("atomic_structure").find("cell")
            for i in range(3):
                A[:, i] = list(map(float, child.find("a{}".format(i + 1)).text.split()))
            A *= Bohr2AA

            species = []
            positions = []
            child = root.find("input").find("atomic_structure").find("atomic_positions")
            for itr in child.iter("atom"):
                species.append(itr.attrib["name"])
                pos = list(map(lambda x: float(x) * Bohr2AA, itr.text.split()))
                positions.append(pos)

            structure = Structure(A, species, positions, coords_are_cartesian=True)
            child = root.find("output").find("total_energy").find("etot")
            ene = hartree2eV * float(child.text)

            Phys = namedtuple("PhysValues", ("energy", "structure"))
            return Phys(np.float64(ene), structure)

    def solver_run_schemes(self):
        return ("mpi_spawn",)
