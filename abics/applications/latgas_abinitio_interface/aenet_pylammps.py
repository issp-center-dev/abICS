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
energy calculator using aenet python interface
"""

from __future__ import annotations

import os.path
from collections import namedtuple
import numpy as np
from pymatgen.core import Structure

from .base_solver import SolverBase, register_solver
from .params import ALParams, DFTParams


def unit_vec(a):
    return a / np.linalg.norm(a)


class AenetPyLammpsSolver(SolverBase):
    """
    Aenetpy solver

    Attributes
    ----------
    path_to_solver : str
        Path to the solver
    input : AenetPySolver.Input
        Input manager
    output : AenetPySolver.Output
        Output manager
    """

    def __init__(self, ignore_species):
        """
        Initialize the solver.

        """
        from mpi4py import MPI
        from lammps import lammps

        super(AenetPyLammpsSolver, self).__init__("")
        self.path_to_solver = self.calculate_energy
        self.input = AenetPyLammpsSolver.Input(ignore_species)
        self.output = AenetPyLammpsSolver.Output()
        self.calc = lammps(
            cmdargs=["-log", "none", "-screen", "none", "-nocite"], comm=MPI.COMM_SELF
        )
        self.rank = MPI.COMM_WORLD.Get_rank()
        
        # debug             cmdargs=["-log", "test{}.log".format(MPI.COMM_WORLD.Get_rank()), "-screen", "none", "-nocite"], comm=MPI.COMM_SELF
    def name(self):
        return "aenetPyLammps"

    def calculate_energy(self, fi, output_dir):
        st = self.input.st
        spec_dict = {}
        for i, sp in enumerate(list(st.symbol_set)):
            spec_dict[sp] = i + 1
        natoms = len(st)
        nspec = len(st.symbol_set)

        latt = st.lattice.matrix
        ax = np.linalg.norm(latt[0])
        bx = np.dot(latt[1], unit_vec(latt[0]))
        by = np.linalg.norm(np.cross(unit_vec(latt[0]), latt[1]))
        cx = np.dot(latt[2], unit_vec(latt[0]))
        cy = (np.dot(latt[1], latt[2]) - bx * cx) / by
        cz = np.sqrt(np.linalg.norm(latt[2]) ** 2.0 - cx ** 2.0 - cy ** 2.0)

        # lmp = lammps(cmdargs=["-log", "none","-nocite"])
        self.calc.command("atom_style atomic")
        self.calc.command("units metal")
        self.calc.command("boundary p p p")

        self.calc.command(f"region box prism 0 {ax} 0 {by} 0 {cz} {bx} {cx} {cy}")
        self.calc.command(f"create_box {nspec} box")
        types = list(map(lambda x: spec_dict[x.name], st.species))
        new_coords = np.dot(st.frac_coords, np.array([[ax,0,0], [bx,by,0], [cx,cy,cz]]))
        self.calc.create_atoms(natoms, None, types, new_coords.ravel())
        
        for i in range(1, nspec + 1):
            self.calc.command(f"mass {i} 1")
        self.calc.commands_list(self.input.pair_pot)
        self.calc.command("run 0")
        ene = self.calc.get_thermo("pe")
        # print(ene)
        self.output.st = st
        self.output.ene = ene
        self.calc.command("clear")

    class Input(object):
        """
        Input manager for Mock

        Attributes
        ----------
        st : pymatgen.Structure
            structure
        """

        st: Structure

        def __init__(self, ignore_species=None):
            self.ignore_species = ignore_species
            # self.st = Structure()

        def from_directory(self, base_input_dir):
            """

            Parameters
            ----------
            base_input_dir : str
                Path to the directory including base input files.
            """
            self.base_input_dir = base_input_dir
            self.pair_pot: list[str] = []
            with open("{}/in.lammps".format(base_input_dir)) as f:
                for line in f:
                    line = line.strip()
                    if line == "" or line.startswith("#"):
                        continue
                    self.pair_pot.append(line)

        def update_info_by_structure(self, structure):
            """
            Update information by structure file

            Parameters
            ----------
            structure : pymatgen.Structure
                Atomic structure
            """
            self.st = structure.copy()
            if self.ignore_species is not None:
                self.st.remove_species(self.ignore_species)

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
            if not os.path.exists(output_dir):
                import shutil

                shutil.copytree(self.base_input_dir, output_dir)

            # self.st.to("POSCAR", os.path.join(output_dir, "structure.vasp"))

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
            Phys = namedtuple("PhysVaules", ("energy", "structure"))
            return Phys(self.ene, self.st)

    def solver_run_schemes(self):
        return ("function",)

    @classmethod
    def create(cls, params: ALParams | DFTParams):
        ignore_species = params.ignore_species
        return cls(ignore_species)


register_solver("aenetpylammps", AenetPyLammpsSolver)
