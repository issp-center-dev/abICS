# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# abICS wrapper of NequIP solver
# Munehiro Kobayashi, Yusuke Konishi (Academeia Co., Ltd.) 2024
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
from typing import Sequence, Dict

import numpy as np
import os, pathlib, shutil, subprocess, shlex
import time

from pymatgen.core import Structure

from abics.util import expand_cmd_path
from abics.applications.latgas_abinitio_interface.base_trainer import TrainerBase
from abics.applications.latgas_abinitio_interface.util import structure_to_XSF

import ase
from ase import io
from ase.calculators.singlepoint import SinglePointCalculator

from nequip.utils import Config
from nequip.scripts import deploy as nequip_deploy


def xsf_to_ase(xsf):
    ase_xsf = ase.io.read(xsf)
    with open(xsf) as f:
        lines = f.readlines()

    tot_energy = float(lines[0].split()[4])
    ase_xsf.calc = SinglePointCalculator(energy=tot_energy, atoms=ase_xsf)
    return ase_xsf

class NequipTrainer(TrainerBase):
    def __init__(
        self,
        structures: Sequence[Structure],
        energies: Sequence[float],
        generate_inputdir: os.PathLike,
        train_inputdir: os.PathLike,
        predict_inputdir: os.PathLike,
        execute_commands: Dict,
        # trainer_type: str,
    ):
        self.structures = structures
        self.energies = energies
        self.generate_inputdir = generate_inputdir
        self.train_inputdir = train_inputdir
        self.predict_inputdir = predict_inputdir
        train_exe = execute_commands["train"]
        self.train_exe = [expand_cmd_path(e) for e in shlex.split(train_exe)]
        self.train_exe.append("input.yaml")
        assert len(self.structures) == len(self.energies)
        self.numdata = len(self.structures)
        self.is_prepared = False
        self.is_trained = False
        self.generate_outputdir = None
        # self.trainer_type = trainer_type

    def prepare(self, latgas_mode = True, st_dir = "nequipXSF"):
        rootdir = os.getcwd()
        xsfdir = os.path.join(rootdir, st_dir)
        
        # prepare XSF files for nequip
        os.makedirs(xsfdir, exist_ok=True)
        os.chdir(xsfdir)
        xsfdir = os.getcwd()
        if latgas_mode:
            for i, st in enumerate(self.structures):
                xsf_string = structure_to_XSF(st, write_force_zero=False)
                xsf_string = (
                    "# total energy = {} eV\n\n".format(self.energies[i]) + xsf_string
                )
                with open("structure.{}.xsf".format(i), "w") as fi:
                    fi.write(xsf_string)
        else:
            for i, st in enumerate(self.structures):
                xsf_string = structure_to_XSF(st, write_force_zero=False)
                xsf_string = (
                    "# total energy = {} eV\n\n".format(self.energies[i]) + xsf_string
                )
                with open("structure.{}.xsf".format(i), "w") as fi:
                    fi.write(xsf_string)

        os.chdir(rootdir)

    def generate_run(self, xsfdir="nequipXSF", generate_dir="generate"):
        # prepare generate
        xsfdir = str(pathlib.Path(xsfdir).resolve())
        if os.path.exists(generate_dir):
            shutil.rmtree(generate_dir)
        # shutil.copytree(self.generate_inputdir, generate_dir)
        os.makedirs(generate_dir, exist_ok=True)
        self.generate_dir = generate_dir
        os.chdir(generate_dir)
        xsf_paths = [
            os.path.join(xsfdir, "structure.{}.xsf".format(i))
            for i in range(self.numdata)
        ]
        ases = [xsf_to_ase(xsf) for xsf in xsf_paths]
        #generate structure.xyz
        ase.io.write("structure.xyz", ases)
        self.generate_outputdir = os.getcwd()
        os.chdir(pathlib.Path(os.getcwd()).parent)
        
    def generate_wait(self):
        interval = 0.1 # sec
        self.is_prepared = False
        if os.path.exists(os.path.join(self.generate_outputdir, "structure.xyz")):
            self.is_prepared = True
        time.sleep(interval)
        if not self.is_prepared:
            raise RuntimeError(f"{self.generate_outputdir}")

    def train(self, train_dir = "train"):
        if not self.is_prepared:
            raise RuntimeError("you have to prepare the trainer before training!")
        if os.path.exists(train_dir):
            shutil.rmtree(train_dir)
        shutil.copytree(self.train_inputdir, train_dir)
        os.chdir(train_dir)

        os.rename(
            os.path.join(self.generate_outputdir, "structure.xyz"),
            os.path.join(os.getcwd(), "structure.xyz"),
        )

        with open(os.path.join(os.getcwd(), "stdout"), "w") as fi:
            subprocess.run(
                self.train_exe, stdout=fi, stderr=subprocess.STDOUT, check=True
            )
        os.chdir(pathlib.Path(os.getcwd()).parent)
        self.is_trained = True

    def new_baseinput(self, baseinput_dir, train_dir = "train"):
        try:
            assert self.is_trained
        except AssertionError as e:
            e.args += "you have to train before getting results!"
    
        baseinput = str(pathlib.Path(baseinput_dir).resolve())
        os.makedirs(baseinput, exist_ok=True)
        shutil.copy(os.path.join(train_dir,"input.yaml"),baseinput)
        os.chdir(train_dir)
        yaml_dic = Config.from_file("input.yaml")
        root = yaml_dic["root"]
        runname = yaml_dic["run_name"]
        nequip_deploy_args = ["build","--train-dir",os.path.join(root,runname),os.path.join(baseinput,"deployed.pth")]
        nequip_deploy.main(nequip_deploy_args)
        os.chdir(pathlib.Path(os.getcwd()).parent)
