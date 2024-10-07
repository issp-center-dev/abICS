from __future__ import annotations

import os
import pathlib
import shlex
import shutil
import subprocess
import time
from typing import Sequence

import ase
from ase.calculators.singlepoint import SinglePointCalculator
from nequip.utils import Config
from pymatgen.core import Structure

from ...util import expand_cmd_path
from . import mlip_3

def to_XSF(structure: Structure, write_force_zero=False):
    """
    Returns a string with the structure in XSF format
    See http://www.xcrysden.org/doc/XSF.html
    """
    lines = []
    app = lines.append

    app("CRYSTAL")
    app("# Primitive lattice vectors in Angstrom")
    app("PRIMVEC")
    cell = structure.lattice.matrix
    for i in range(3):
        app(" %.14f %.14f %.14f" % tuple(cell[i]))

    cart_coords = structure.cart_coords
    app("# Cartesian coordinates in Angstrom.")
    app("PRIMCOORD")
    app(" %d 1" % len(cart_coords))
    species = structure.species
    site_properties = structure.site_properties
    if "forces" not in site_properties.keys():
        write_force_zero = True
    else:
        forces = site_properties["forces"]

    if write_force_zero:
        for a in range(len(cart_coords)):
            app(
                str(species[a])
                + " %20.14f %20.14f %20.14f" % tuple(cart_coords[a])
                + " 0.0 0.0 0.0"
            )
    else:
        for a in range(len(cart_coords)):
            app(
                str(species[a])
                + " %20.14f %20.14f %20.14f" % tuple(cart_coords[a])
                + " %20.14f %20.14f %20.14f" % tuple(forces[a])
            )

    return "\n".join(lines)

class mlip_3_trainer:
    def __init__(
        self,
        structures: Sequence[Structure],
        energies: Sequence[float],
        generate_inputdir: os.PathLike,
        train_inputdir: os.PathLike,
        predict_inputdir: os.PathLike,
        generate_exe: str,
        train_exe: str,
    ):
        self.structures = structures
        self.energies = energies
        self.generate_inputdir = generate_inputdir
        self.train_inputdir = train_inputdir
        self.predict_inputdir = predict_inputdir
        self.train_exe = [
            expand_cmd_path(e) for e in shlex.split(train_exe)
        ]
        self.train_exe += ["train", "input.almtp", "input.cfg", "--save_to=./pot.almtp", "--iteration_limit=100", "--al_mode=nbh"]
        assert len(self.structures) == len(self.energies)
        self.numdata = len(self.structures)
        self.is_prepared = False
        self.is_trained = False
        self.generate_outputdir = None
        self.latgas_mode = True

    def prepare(self, latgas_mode=True, st_dir="mlip_3_XSF"):
        rootdir = os.getcwd()
        xsfdir = os.path.join(rootdir, st_dir)

        # prepare XSF files
        os.makedirs(xsfdir, exist_ok=True)
        os.chdir(xsfdir)
        xsfdir = os.getcwd()
        if latgas_mode:
            for i, st in enumerate(self.structures):
                xsf_string = to_XSF(st, write_force_zero=False)
                xsf_string =\
                    f"# total energy = {self.energies[i]} eV\n\n{xsf_string}"
                with open(f"structure.{i}.xsf", "w") as fi:
                    fi.write(xsf_string)
        else:
            for i, st in enumerate(self.structures):
                xsf_string = to_XSF(st, write_force_zero=False)
                xsf_string =\
                    f"# total energy = {self.energies[i]} eV\n\n{xsf_string}"
                with open(f"structure.{i}.xsf", "w") as fi:
                    fi.write(xsf_string)

        self.latgas_mode = latgas_mode
        os.chdir(rootdir)

    def generate_run(self, xsfdir="mlip_3_XSF", generate_dir="generate"):
        # prepare generate
        cfgdir = str(pathlib.Path(xsfdir).resolve())
        if os.path.exists(generate_dir):
            shutil.rmtree(generate_dir)
        shutil.copytree(xsfdir, generate_dir)
        os.chdir(generate_dir)

        # prepare CFG file for MLIP-3
        if self.latgas_mode:
            cfg_string = ""
            for i, st in enumerate(self.structures):
                lines = mlip_3.to_CFG(st, self.energies[i], write_force_zero=False)
                cfg_string = cfg_string + lines
            with open(f"input.cfg", "w") as fi:
                fi.write(cfg_string)           
        else:
            cfg_string = ""
            for i, st in enumerate(self.structures):
                lines = mlip_3.to_CFG(st, self.energies[i], write_force_zero=False)
                cfg_string = cfg_string + lines
            with open(f"input.cfg", "w") as fi:
                fi.write(cfg_string)           

        self.generate_outputdir = os.getcwd()
        os.chdir(pathlib.Path(os.getcwd()).parent)

    def generate_wait(self):
        interval = 0.1  # sec
        #self.is_prepared = False
        #if os.path.exists(
        #    os.path.join(self.generate_outputdir, "input.cfg")
        #):
        #    self.is_prepared = True
        self.is_prepared = True
        time.sleep(interval)
        if not self.is_prepared:
            raise RuntimeError(f"{self.generate_outputdir}")

    def train(self, train_dir="train"):
        if not self.is_prepared:
            raise RuntimeError(
                "you have to prepare the trainer before training!"
            )
        if os.path.exists(train_dir):
            shutil.rmtree(train_dir)
        shutil.copytree(self.generate_outputdir, train_dir)
        shutil.copy(os.path.join(self.train_inputdir, "input.almtp"), train_dir)
        os.chdir(train_dir)
        #command = self.train_exe + " train input.almtp input.cfg --save_to=./out/pot.mtp --interaction_limit=100 --al_mode=nbh"
        command = self.train_exe
        #print(os.getcwd())
        #print(command)
        #print(os.path.exists("input.cfg"))

        with open(os.path.join(os.getcwd(), "stdout"), "w") as fi:
            subprocess.run(
                self.train_exe, stdout=fi, stderr=subprocess.STDOUT, check=True
            )
        os.chdir(pathlib.Path(os.getcwd()).parent)
        self.is_trained = True

    def new_baseinput(self, baseinput_dir, train_dir="train"):
        try:
            assert self.is_trained
        except AssertionError as e:
            e.args += "you have to train before getting results!"

        baseinput = str(pathlib.Path(baseinput_dir).resolve())
        os.makedirs(baseinput, exist_ok=True)
        shutil.copy(os.path.join(train_dir, "input.cfg"), baseinput)
        shutil.copy(os.path.join(train_dir, "pot.almtp"), baseinput)

