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

from .base_solver import SolverBase
from collections import namedtuple
import numpy as np
import os
import scipy.constants as spc
import subprocess
from pymatgen import Structure

hartree2eV = spc.value("Hartree energy in eV")
Bohr2AA = spc.value("Bohr radius") * 1e10


class OpenMXInputFile(dict):
    """ OpenMX Input

    This is a dictionary storing OpenMX input parameters
    with additionary information, vec_list

    Attributes
    ----------
    vec_list: list
        vectors
    """

    def __init__(self, input_file):
        """

        Parameters
        ----------
        input_file: str
            inputfile
        """
        super().__init__(self)
        self.vec_list = [
            "Atoms.UnitVectors",
            "Atoms.SpeciesAndCoordinates",
            "MD.Fixed.XYZ",
        ]
        with open(input_file, "r") as f:
            lines = f.readlines()
            # delete comment out
            list_flag = False
            vals = []
            for line in lines:
                line = line.strip().split("#")[0]
                if len(line) >= 2 or line != "":
                    words = line.split()
                    if words[0][0] == "<":
                        list_flag = True
                        vals = []
                    elif words[0][-1:] == ">":
                        self[line[:-1]] = vals
                        self.vec_list.append(line[:-1])
                        list_flag = False
                    else:
                        if list_flag:
                            vals.append(words)
                        else:
                            self[words[0]] = words[1:]
                    self.vec_list = list(set(self.vec_list))


class OpenMXSolver(SolverBase):
    """
    This class defines the OpenMX solver.
    """

    def __init__(self, path_to_solver):
        """
        Initialize the solver.

        Parameters
        ----------
        path_to_solver : str
                      Path to the solver.
        """
        super(OpenMXSolver, self).__init__(path_to_solver)
        self.path_to_solver = path_to_solver
        self.input = OpenMXSolver.Input(self.path_to_solver)
        self.output = OpenMXSolver.Output(self.input)

    def name(self):
        """
        Return the solver name.

        Returns
        -------
        "OpenMX"

        """
        return "OpenMX"

    class Input(object):
        def __init__(self, path_to_solver):
            self.base_info = None
            self.pos_info = None
            self.base_openmx_input = None
            self.path_to_solver = path_to_solver

        def cleanup(self, rundir):
            """
            Clean up the directory where the solver runs.

            Parameters
            ----------
            rundir: str
                Name of the directory where the solver runs.
            -------

            """

            checkfile = os.path.join(rundir, self.datadir, self.filetocheck)
            if os.path.exists(checkfile):
                os.remove(checkfile)

        def from_directory(self, base_input_dir):
            """

            Get base information for OpenMX from the input directory.

            Parameters
            ----------
            base_input_dir: str
                Name of the directory where the base input file is located.

            Returns
            -------
            self.base_openmx_input:  OpenMXInput
            """
            # TODO
            # check the base input file name (now, set "base.dat")
            self.base_openmx_input = OpenMXInputFile(
                os.path.join(os.getcwd(), base_input_dir, "base.dat")
            )
            self.vps_info = self._get_vps_info(self.base_openmx_input)

            return self.base_openmx_input

        def _get_vps_info(self, openmx_input):
            """
            Get vps information from the DFT_DATA19 directory.

            Parameters
            ----------
            openmx_input: dict
                Check whether Data.path is included in the base.dat file

            Returns
            -------
            vps_info:  dict (key: specie, value: electron number)
                Dictionary for vps_info.

            """
            vps_info = {}
            if "DATA.PATH" in openmx_input:
                path = openmx_input["DATA.PATH"][0]
            else:
                cmd = "which {}".format(self.path_to_solver)
                fullpath = subprocess.check_output(cmd.split()).splitlines()[0].decode()
                path = os.path.join(os.path.dirname(fullpath), "..", "DFT_DATA19")
                openmx_input["DATA.PATH"] = [path]
            with open(os.path.join(path, "vps_info.txt"), "r") as f:
                lines = f.readlines()
                for line in lines:
                    vps_info[line.split()[0]] = float(line.split()[1])
            return vps_info

        def update_info_by_structure(self, structure):
            """
            Update base_openmx_input by structure

            Parameters
            ----------
            structure: pymatgen.core.Structure
                Structure for getting atom's species and coordinates

            """
            # Get lattice information
            A = structure.lattice.matrix
            # Update unitvector information
            self.base_openmx_input["Atoms.UnitVectors.Unit"] = ["Ang"]
            self.base_openmx_input["Atoms.UnitVectors"] = A  # numpy.ndarray
            nat = len(structure.sites)
            # For write_input, Atoms.Number must be list format.
            self.base_openmx_input["Atoms.Number"] = [nat]
            self.base_openmx_input["Atoms.SpeciesAndCoordinates.Unit"] = ["FRAC"]
            self.base_openmx_input["Atoms.SpeciesAndCoordinates"] = [
                [0, "", 0.0, 0.0, 0.0, 0.0, 0.0]
            ] * nat

            mag = [0.0] * nat
            if "magnetization" in structure.site_properties:
                mag = structure.site_properties["magnetization"]

            # Get electron_info
            atomic_species = self.base_openmx_input["Definition.of.Atomic.Species"]
            for idx, site in enumerate(structure.sites):
                electron_number = self.vps_info[
                    [
                        specie[2]
                        for specie in atomic_species
                        if specie[0] == str(site.specie)
                    ][0]
                ]
                self.base_openmx_input["Atoms.SpeciesAndCoordinates"][idx] = [
                    idx + 1,
                    site.specie,
                    site.a,
                    site.b,
                    site.c,
                    0.5 * electron_number + mag[idx],
                    0.5 * electron_number - mag[idx],
                ]
            if "seldyn" in structure.site_properties:
                seldyn_arr = structure.site_properties["seldyn"]
                self.base_openmx_input["MD.Fixed.XYZ"] = [
                    [0, int(False), int(False), int(False)]
                ] * nat
                for idx, dyn_info in enumerate(seldyn_arr):
                    fix_info = (~np.array(dyn_info)).astype(int)
                    self.base_openmx_input["MD.Fixed.XYZ"][idx] = [
                        idx + 1,
                        fix_info[0],
                        fix_info[1],
                        fix_info[2],
                    ]

        def write_input(self, output_dir):
            """
            Generate input files at each directory.

            Parameters
            ----------
            output_dir: str
                Path to the output directory
            """
            try:
                os.makedirs(output_dir)
            except:
                pass
            output_file = os.path.join(
                output_dir, "{}.dat".format(self.base_openmx_input["System.Name"][0])
            )
            self.base_openmx_input['System.CurrrentDirectory'] = [output_dir + '/']
            with open(output_file, "w") as f:
                for key, values in self.base_openmx_input.items():
                    if key in self.base_openmx_input.vec_list:
                        print_stamp = "<{}\n".format(key)
                        for value_list in values:
                            for value in value_list:
                                print_stamp += " {}".format(value)
                            print_stamp += "\n"
                        print_stamp += "{}>\n".format(key)
                    else:
                        print_stamp = key
                        for value_list in values:
                            print_stamp += " {}".format(value_list)
                        print_stamp += "\n"
                    f.write(print_stamp)

        def cl_args(self, nprocs, nthreads, output_dir):
            """
            Make command line to execute OpenMX from args.

            Parameters
            ----------
            nprocs: int
                Number of processes (not used).

            nthreads: int
                Number of threads.

            output_dir: str
                Output directory.

            Returns
            -------
            clargs: dict
                command line arguments
            """
            clargs = [
                "{}.dat".format(
                    os.path.join(output_dir, self.base_openmx_input["System.Name"][0])
                ),
                "-nt",
                str(nthreads),
            ]
            return clargs

        # def submit: Use submit defined in run_base_mpi.py

    class Output(object):
        def __init__(self, input):
            self.input = input
            pass

        def get_results(self, output_dir):
            """
            Read results from files in output_dir and calculate values

            Parameters
            ----------
            output_dir: str
                    Output directory name

            Returns
            -------
            Phys: namedtuple
                results

            """
            # Read results from files in output_dir and calculate values
            output_file = os.path.join(
                output_dir,
                "{}.out".format(self.input.base_openmx_input["System.Name"][0]),
            )
            with open(output_file, "r") as fi:
                lines = fi.readlines()
                lines_strip = [line.strip() for line in lines]
                # Get total energy
                Utot = float(
                    [line for line in lines_strip if "Utot." in line][0].split()[1]
                )
                # Change energy unit from Hartree to eV
                Utot *= hartree2eV

            # Get Cell information from dat# file
            A = np.zeros((3, 3))
            output_file = os.path.join(
                output_dir,
                "{}.dat#".format(self.input.base_openmx_input["System.Name"][0]),
            )
            with open(output_file, "r") as fi:
                lines = fi.readlines()
                lines_strip = [line.strip() for line in lines]
                # Get Cell information
                # Read Atoms.UnitVectors.Unit
                Atoms_UnitVectors_Unit = [
                    line.split()[1]
                    for line in lines_strip
                    if "Atoms.UnitVectors.Unit" in line
                ][0]
                line_number_unit_vector_start = [
                    i
                    for i, line in enumerate(lines_strip)
                    if "<Atoms.UnitVectors" in line
                ][0]
                for i, line in enumerate(
                    lines_strip[
                        line_number_unit_vector_start
                        + 1 : line_number_unit_vector_start
                        + 4
                    ]
                ):
                    A[:, i] = list(line.split())
                if Atoms_UnitVectors_Unit == "AU":
                    A *= Bohr2AA

            # Note:
            # Since xyz format of OpenMX is not correct (3 columns are added at each line).
            # pymatgen.io.xyz can not work.
            output_file = os.path.join(
                output_dir,
                "{}.xyz".format(self.input.base_openmx_input["System.Name"][0]),
            )
            species = []
            positions = []
            with open(output_file, "r") as fi:
                # Get atomic species and positions
                lines = fi.readlines()
                lines_strip = [line.strip().split() for line in lines]
                for line in lines_strip[2:]:
                    species.append(line[0])
                    positions.append([float(s) for s in line[1:4]])
            structure = Structure(A, species, positions, coords_are_cartesian=True)
            Phys = namedtuple("PhysValues", ("energy", "structure"))
            return Phys(np.float64(Utot), structure)

    def solver_run_schemes(self):
        return ("mpi_spawn_ready", "mpi_spawn_wrapper")
