# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2020- The University of Tokyo
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
from typing import MutableMapping, Any

import sys, datetime

import os, sys
import threading
import itertools

import numpy as np
import networkx as nx
from pymatgen.core import Structure

from abics import __version__
from abics.applications.latgas_abinitio_interface.params import DFTParams, TrainerParams
from abics.applications.latgas_abinitio_interface import aenet_trainer
from abics.applications.latgas_abinitio_interface import map2perflat
from abics.applications.latgas_abinitio_interface.naive_matcher import naive_mapping
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)

import logging
import extra.loggers
from pathlib import Path

extra.loggers.set_log_handles(app_name="train", level=logging.INFO, mpi_log="master")
logger = logging.getLogger("main")


def main_impl(params_root: MutableMapping):
    if not os.path.exists("ALloop.progress"):
        # print("abics_mlref has not run yet.")
        logger.error("abics_mlref has not run yet.")
        sys.exit(1)
    ALdirs = []
    with open("ALloop.progress", "r") as fi:
        lines = fi.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith("AL"):
                ALdirs.append(line)
        last_li = lines[-1].strip()
    if not last_li.startswith("AL"):
        # print("You shouldn't run train now.")
        # print("Either abics_sampling or abics_mlref first.")
        logger.error("You shouldn't run train now.")
        logger.error("Either abics_sampling or abics_mlref first.")
        sys.exit(1)

    dftparams = DFTParams.from_dict(params_root["sampling"]["solver"])
    ensemble = dftparams.ensemble
    base_input_dir = dftparams.base_input_dir

    trainerparams = TrainerParams.from_dict(params_root["train"])
    ignore_species = trainerparams.ignore_species
    trainer_commands = trainerparams.exe_command
    trainer_type = trainerparams.solver
    trainer_input_dirs = trainerparams.base_input_dir

    configparams = DFTConfigParams.from_dict(params_root["config"])
    config = defect_config(configparams)
    species = config.structure.symbol_set
    dummy_sts = {sp: config.dummy_structure_sp(sp) for sp in species}

    if trainer_type != "aenet":
        # print("Unknown trainer: ", trainer_type)
        logger.error("Unknown trainer: ", trainer_type)
        sys.exit(1)

    rootdir = os.getcwd()
    structures = []
    energies = []


    # print("-Mapping relaxed structures in AL* to on-lattice model...")
    # sys.stdout.flush()
    logger.info("-Mapping relaxed structures in AL* to on-lattice model...")
    # val_map is a list of list [[sp0, vac0], [sp1, vac1], ...]
    # if vac_map:
    #    vac_map = {specie: vacancy for specie, vacancy in vac_map}
    # else:
    #    vac_map = {}

    # we first group species that share sublattices together

    G = nx.Graph()
    G.add_nodes_from(species)
    for sublattice in config.defect_sublattices:
        groups = sublattice.groups
        sp_list = []
        for group in groups:
            sp_list.extend(group.species)
        for pair in itertools.combinations(sp_list, 2):
            G.add_edge(*pair)
    sp_groups = nx.connected_components(G)
    # print(list(sp_groups))
    # sys.exit(0)
    dummy_sts_share : list[tuple[Structure, list]] = []
    for c in nx.connected_components(G):
        # merge dummy structures for species that share sublattices
        sps = list(c)

        coords = np.concatenate([dummy_sts[sp].frac_coords for sp in sps], axis=0)
        st_tmp = Structure(
            dummy_sts[sps[0]].lattice,
            species=["X"] * coords.shape[0],
            coords=coords,
        )
        st_tmp.merge_sites(mode="delete")
        dummy_sts_share.append((st_tmp, sps))

    # for i,st in enumerate(dummy_sts_share):
    #    st.to(fmt="POSCAR", filename="{}.dummy.vasp".format(i))

    for dir in ALdirs:
        rpl = 0
        while os.path.isdir(os.path.join(dir, str(rpl))):
            os.chdir(os.path.join(dir, str(rpl)))
            energies_ref = []
            step_ids = []
            with open("energy_corr.dat") as fi:
                for line in fi:
                    words = line.split()
                    energies_ref.append(float(words[1]))
                    step_ids.append(int(words[2]))
            for step_id, energy in zip(step_ids, energies_ref):
                structure: Structure = Structure.from_file(f"structure.{step_id}.vasp")
                mapped_sts = []
                mapping_success = True
                for dummy_st, specs in dummy_sts_share:
                    # perform sublattice by sublattice mapping
                    sp_rm = list(filter(lambda s: s not in specs, species))
                    st_tmp = structure.copy()
                    st_tmp.remove_species(sp_rm)
                    num_sp = len(st_tmp)
                    # map to perfect lattice for this species
                    st_tmp = map2perflat(dummy_st, st_tmp)
                    st_tmp.remove_species(["X"])
                    mapped_sts.append(st_tmp)
                    if num_sp != len(st_tmp):
                        # print(
                        #     f"--mapping failed for structure {step_id} in replica {rpl}"
                        # )
                        logger.info(f"--mapping failed for structure {step_id} in replica {rpl}")
                        mapping_success = False

                for sts in mapped_sts[1:]:
                    for i in range(len(sts)):
                        mapped_sts[0].append(sts[i].species_string, sts[i].frac_coords)
                if ignore_species:
                    mapped_sts[0].remove_species(ignore_species)
                if mapping_success:
                    structures.append(mapped_sts[0])
                    energies.append(energy)
            rpl += 1
            os.chdir(rootdir)
    # print("--Finished mapping")
    # sys.stdout.flush()
    logger.info("--Finished mapping")
    
    generate_input_dirs = []
    train_input_dirs = []
    predict_input_dirs = []

    if dftparams.ensemble:
        if len(trainer_input_dirs) != len(base_input_dir):
            # print(
            #     "You must set the number of trainer input dirs equal to baseinput dirs for ensemble NNP"
            # )
            logger.error("You must set the number of trainer input dirs equal to baseinput dirs for ensemble NNP")
            sys.exit(1)
    for d in trainer_input_dirs:
        generate_input_dirs.append(os.path.join(d, "generate"))
        train_input_dirs.append(os.path.join(d, "train"))
        predict_input_dirs.append(os.path.join(d, "predict"))

    generate_exe = trainer_commands[0]
    train_exe = trainer_commands[1]

    trainers = []
    for i in range(len(trainer_input_dirs)):
        trainers.append(
            aenet_trainer(
                structures,
                energies,
                generate_input_dirs[i],
                train_input_dirs[i],
                predict_input_dirs[i],
                generate_exe,
                train_exe,
            )
        )

    trainers[0].prepare()
    # We use threads to parallelize generate.x over ensemble members
    """threads = []
    for i, trainer in enumerate(trainers):
        threads.append(
            threading.Thread(
                target=trainer.generate(generate_dir="generate{}".format(i))
                )
            )
        threads[-1].start()
    for t in threads:
        t.join()
    """
    for i, trainer in enumerate(trainers):
        # print(f"-Running generate run in generate{i}")
        # sys.stdout.flush()
        logger.info(f"-Running generate run in generate{i}")
        trainer.generate_run(generate_dir="generate{}".format(i))

    for trainer in trainers:
        trainer.generate_wait()
    # print(f"--Finished generate run(s)")
    logger.info(f"--Finished generate run(s)")
    
    # We use MPI version of train.x so no need to write parallel code here
    for i, trainer in enumerate(trainers):
        # print(f"-Training run in train{i}")
        # sys.stdout.flush()
        logger.info(f"-Training run in train{i}")
        trainer.train(train_dir="train{}".format(i))
        # print(f"--Training run finished in train{i}")
        # print(f"-Preparing NN model for abics_sampling in {base_input_dir[i]}")
        # sys.stdout.flush()
        logger.info(f"--Training run finished in train{i}")
        logger.info(f"-Preparing NN model for abics_sampling in {base_input_dir[i]}")
        trainer.new_baseinput(base_input_dir[i])
        # print(f"--Success.")
        logger.info(f"--Success.")

    with open("ALloop.progress", "a") as fi:
        # print("-Writing ALloop.progress")
        logger.info("-Writing ALloop.progress")
        fi.write("train\n")
        fi.flush()
        os.fsync(fi.fileno())
    now = datetime.datetime.now()
    # print("-Let's run abics_sampling next!")
    # print(f"Exiting normally on {now}.\n")
    logger.info("-Let's run abics_sampling next!")
    logger.info(f"Exiting normally on {now}.\n")


def main():
    import toml
    now = datetime.datetime.now()
    # print(f"Running abics_train (abICS v{__version__}) on {now}")
    logger.info(f"Running abics_train (abICS v{__version__}) on {now}")
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    # print(f"-Reading input from: {tomlfile}")
    logger.info(f"-Reading input from: {tomlfile}")
    params_root = toml.load(tomlfile)
    main_impl(params_root)


if __name__ == "__main__":
    main()
