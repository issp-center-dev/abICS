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
To deal with QuantumESPRESSO

"""

from collections import namedtuple
import xml.etree.ElementTree as ET
import operator
import os
import os.path

import numpy as np
import scipy.constants as spc
from pymatgen import Structure
from qe_tools.parsers import PwInputFile

from .base_solver import SolverBase
from ...util import expand_path
from ...exception import InputError


hartree2eV = spc.value("Hartree energy in eV")
Bohr2AA = spc.value("Bohr radius") * 1e10


class QESolver(SolverBase):
    """
    Solver class dealing with QuantumESPRESSO (with new XML format).
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
        """
        Returns
        -------
        solver_name : str
            Name of solver, "QuantumESPRESSO".
        """
        return "QuantumESPRESSO"

    class Input(object):
        """
        QE input files manager.

        Attributes
        ----------
        pwi : qe_tools.parsers.PwInputFile
            QE input.
        datadir : str
            Path to the data directory of QE.
        filetocheck : str
            Name of the file to be used for check finished.
        """

        def __init__(self):
            self.pwi = None
            self.datadir = "pwscf.save"
            self.filetocheck = "data-file-schema.xml"

        def cleanup(self, workdir):
            """
            Remove file for check finished.

            Parameters
            ----------
            workdir : str
                Path to the working directory.
            """
            checkfile = os.path.join(workdir, self.datadir, self.filetocheck)
            if os.path.exists(checkfile):
                os.remove(checkfile)

        def check_finished(self, workdir):
            """
            Check if the solver program has finished or not.

            Parameters
            ----------
            workdir : str
                Path to the working directory.
            """
            f = os.path.join(workdir, self.datadir, self.filetocheck)
            if not os.path.exists(f):
                return False
            try:
                tree = ET.parse(f)
            except ET.ParseError:
                return False

            root = tree.getroot()
            control = root.find("input").find("control_variables")
            calculation = control.find("calculation").text
            if calculation == "scf":
                return True
            elif calculation == "relax":
                scf_conv = root.find("output").find("convergence_info").find("scf_conv")
                scf_converged = scf_conv.find("convergence_achieved").text.lower() == "true"
                # scf_maxstep = int(root.find("input").find("electron_control").find("max_nstep").text)
                # scf_nstep = int(scf_conv.find("n_scf_steps").text)

                if not scf_converged:
                    # scf does not converged and QE stopped
                    return True

                opt_maxstep = int(control.find("nstep").text)
                opt_conv = root.find("output").find("convergence_info").find("opt_conv")
                opt_converged = opt_conv.find("convergence_achieved").text.lower() == "true"
                if opt_converged:
                    return True
                opt_nstep = int(opt_conv.find("n_opt_steps").text)
                return opt_nstep == opt_maxstep
            else:
                raise InputError("calculation is {}, but this is not yet supported".format(calculation))

        def from_directory(self, base_input_dir):
            """
            Initialize information from files in base_input_dir.

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            """

            inputfile = os.path.join(os.getcwd(), base_input_dir, "scf.in")
            self.pwi = PwInputFile(inputfile)
            for name in ['CONTROL', 'SYSTEM']:
                if name not in self.pwi.namelists:
                    raise InputError("{} cannot found in {}".format(name, inputfile))
            calculation = self.pwi.namelists["CONTROL"]["calculation"]
            if calculation not in ("scf", "relax"):
                raise InputError("Sorry, {} is not yet supported".format(calculation))
            for name in ['ELECTRONS', 'IONS', 'CELL']:
                if name not in self.pwi.namelists:
                    self.pwi.namelists[name] = {}
            self.pwi.namelists["CONTROL"]["prefix"] = "pwscf"
            self.pwi.namelists["CONTROL"]["pseudo_dir"] = expand_path(
                self.pwi.namelists["CONTROL"]["pseudo_dir"], os.getcwd()
            )

        def update_info_by_structure(self, structure):
            """
            Update information by atomic structure.

            Parameters
            ----------
            structure : pymatgen.Structure
                Atomic structure
            """
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

            if "seldyn" in structure.site_properties:
                seldyn_arr = structure.site_properties["seldyn"]
                for idx, dyn_info in enumerate(seldyn_arr):
                    self.pwi.atomic_positions["fixed_coords"][idx] = list(
                        map(operator.not_, dyn_info)
                    )

        def update_info_from_files(self, output_dir, rerun):
            """
            Do nothing.
            """
            pass

        def write_input(self, output_dir):
            """
            Generate input files of the solver program.

            Parameters
            ----------
            output_dir : str
                Path to working directory.
            """
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
                if "points" in self.pwi.k_points:
                    for kp in self.pwi.k_points["points"]:
                        f.write(" {}".format(kp))
                    for offset in self.pwi.k_points["offset"]:
                        f.write(" {}".format(int(2 * offset)))
                    f.write("\n")

        def file_for_check_finished(self, workdir):
            return os.path.join(workdir, self.datadir, self.filetocheck)

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
            return ["-in", os.path.join(workdir, "scf.in")]

    class Output(object):
        """
        Output manager.
        """

        def __init__(self, prefix):
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
                and coordinates is measured in the units of Angstrom.
            """
            # Read results from files in output_dir and calculate values
            tree = ET.parse(os.path.join(workdir, "pwscf.save", "data-file-schema.xml"))
            root = tree.getroot()
            A = np.zeros((3, 3))
            child = root.find("input").find("atomic_structure").find("cell")
            for i in range(3):
                A[:, i] = list(map(float, child.find("a{}".format(i + 1)).text.split()))
            A *= Bohr2AA

            species = []
            positions = []
            output = root.find("output")
            child = output.find("atomic_structure").find("atomic_positions")
            for itr in child.iter("atom"):
                species.append(itr.attrib["name"])
                pos = list(map(lambda x: float(x) * Bohr2AA, itr.text.split()))
                positions.append(pos)

            structure = Structure(A, species, positions, coords_are_cartesian=True)
            child = output.find("total_energy").find("etot")
            ene = hartree2eV * float(child.text)

            Phys = namedtuple("PhysValues", ("energy", "structure"))
            return Phys(np.float64(ene), structure)

    def solver_run_schemes(self):
        """
        Returns
        -------
        schemes : tuple[str]
            Implemented runner schemes.
        """
        return "mpi_spawn"

