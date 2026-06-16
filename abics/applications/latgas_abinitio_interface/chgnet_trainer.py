# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# abICS wrapper of CHGNet solver
# Masashi Noda, Yusuke Konishi (Academeia Co., Ltd.) 2025
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
from ase.calculators.morse import MorsePotential

from chgnet.model.model import CHGNet
from chgnet.data.dataset import StructureData, get_train_val_test_loader
from chgnet.trainer import Trainer

import yaml

def xsf_to_ase(xsf):
    ase_xsf = ase.io.read(xsf)
    with open(xsf) as f:
        lines = f.readlines()

    tot_energy = float(lines[0].split()[4])
    ase_xsf.calc = SinglePointCalculator(energy=tot_energy, atoms=ase_xsf)
    ase_xsf.calc = MorsePotential()
    forces = ase_xsf.get_forces()
    ase_xsf.calc = SinglePointCalculator(energy=tot_energy, forces=forces, atoms=ase_xsf)
    return ase_xsf

class CHGNetTrainer(TrainerBase):
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
        #train_exe = execute_commands["train"]
        #self.train_exe = [expand_cmd_path(e) for e in shlex.split(train_exe)]
        assert len(self.structures) == len(self.energies)
        self.numdata = len(self.structures)
        self.is_prepared = False
        self.is_trained = False
        self.generate_outputdir = None
        self.chgnet_params = {
            "finetuning" : True,
            "batch_size" : 4, 
            "train_ratio" : 0.9, 
            "val_ratio" : 0.05,
            "learning_rate" : 1e-2,
            "epochs" : 5,
            "model_params" : {
                "atom_fea_dim" : 64
            } 
        }
        yaml_file = os.path.join(train_inputdir, "input.yaml")
        with open(yaml_file, "r") as f:
            self.chgnet_params.update(yaml.safe_load(f))
        if "mlp_hidden_dims" in self.chgnet_params["model_params"].keys():
            self.chgnet_params["model_params"]["mlp_hidden_dims"] = tuple(self.chgnet_params["model_params"]["mlp_hidden_dims"])
        # self.trainer_type = trainer_type

    def prepare(self, latgas_mode = True, st_dir = "chgnetXSF"):
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

    def generate_run(self, xsfdir="chgnetXSF", generate_dir="generate"):
        # prepare generate
        xsfdir = str(pathlib.Path(xsfdir).resolve())
        if os.path.exists(generate_dir):
            shutil.rmtree(generate_dir)
        shutil.copytree(self.generate_inputdir, generate_dir)
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

        # read structure.xyz as ase
        atoms_list = io.read("structure.xyz", index=":")
        # make Structure list and energy list
        structures = []
        energies = []
        forces = []

        for atoms in atoms_list:
            # Convert ASE Atoms -> Pymatgen Structure
            structure = Structure(
                lattice=atoms.get_cell(),
                species=atoms.get_chemical_symbols(),
                coords=atoms.get_positions(),
                coords_are_cartesian=True
            )
            structures.append(structure)

            comp = structure.composition.as_dict()
            num_atoms = sum([comp[key] for key in comp.keys()])

            # Get energy and force
            energy = atoms.get_potential_energy()
            force = atoms.get_forces()
            energies.append(energy/num_atoms)
            forces.append(force)

        dataset = StructureData(
            structures=structures,
            energies=energies,
            forces = forces,
        )
        train_loader, val_loader, test_loader = get_train_val_test_loader(
            dataset, 
            batch_size=self.chgnet_params["batch_size"], 
            train_ratio=self.chgnet_params["train_ratio"], 
            val_ratio=self.chgnet_params["val_ratio"]
        )
        if self.chgnet_params["finetuning"]:
            chgnet = CHGNet(atom_graph_cutoff=7.5, bond_graph_cutoff=6.0)
            chgnet = chgnet.load()
        else:
            chgnet = CHGNet(
                **self.chgnet_params["model_params"]
            )
        
        trainer = Trainer(
            model=chgnet,
            targets="ef",
            force_loss_ratio=0.0,
            optimizer="Adam",
            criterion="MSE",
            learning_rate=self.chgnet_params["learning_rate"],
            epochs=self.chgnet_params["epochs"],
        )
        trainer.train(train_loader, val_loader, test_loader, save_dir="chgnet_out", train_composition_model=True)

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
        # Search for bestE_*.pth.tar in the chgnet_out directory and copy the 0th as deployed.pth.tar
        bestE_files = [f for f in os.listdir("chgnet_out") if f.startswith("bestE")]
        shutil.copy(os.path.join("chgnet_out",bestE_files[0]),os.path.join(baseinput,"deployed.pth.tar"))
        os.chdir(pathlib.Path(os.getcwd()).parent)
