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

from abics.replica_params import RefParams
from abics.applications.latgas_abinitio_interface.params import DFTParams, TrainerParams
from abics.applications.latgas_abinitio_interface import aenet_trainer
from abics.applications.latgas_abinitio_interface import map2perflat
from abics.applications.latgas_abinitio_interface.naive_matcher import naive_mapping
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)

from pymatgen import Structure
import numpy as np
import os, sys


def main_impl(tomlfile):
    if not os.path.exists("ALloop.progress"):
        print("You shouldn't run train now. Run activelearn first.")
    with open("ALloop.progress", "r") as fi:
        last_li = fi.readlines()[-1]
    if "AL" not in last_li:
        print("You shouldn't run train now. Either MC or activelearn first.")
        sys.exit(1)

    # We need info from MC steps as well as parameters
    # specific to training
    rxparams = RefParams.from_toml(tomlfile)
    nreplicas = rxparams.nreplicas
    nsteps = rxparams.nsteps
    sample_frequency = rxparams.sample_frequency
    sample_ids = list(range(0, nsteps, sample_frequency))
    
    DFTparams = DFTParams.from_toml(tomlfile)
    base_input_dir = DFTparams.base_input_dir

    trainerparams = TrainerParams.from_toml(tomlfile)
    ignore_species = trainerparams.ignore_species
    trainer_commands = trainerparams.exe_command
    trainer_type = trainerparams.solver
    trainer_input_dirs = trainerparams.base_input_dir
    
    configparams = DFTConfigParams.from_toml(tomlfile)
    config = defect_config(configparams)
    species = config.structure.symbol_set
    dummy_sts = {sp: config.dummy_structure_sp(sp) for sp in species}

    if trainer_type != "aenet":
        print("Unknown trainer: ", trainer_type)
        sys.exit(1)

    ls = os.listdir()
    rootdir = os.getcwd()
    ALdirs = [dir for dir in ls if "AL" in dir and os.path.isdir(dir)]
    structures = []
    energies = []

    # val_map is a list of list [[sp0, vac0], [sp1, vac1], ...]
    #if vac_map:
    #    vac_map = {specie: vacancy for specie, vacancy in vac_map}
    #else:
    #    vac_map = {}

    # we first group species that share sublattices together
    import networkx as nx
    import itertools
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
    #print(list(sp_groups))
    #sys.exit(0)
    dummy_sts_share = []
    for c in nx.connected_components(G):
        # merge dummy structures for species that share sublattices
        sps = list(c)

        coords = np.concatenate([dummy_sts[sp].frac_coords for sp in sps], axis = 0)
        st_tmp = Structure(dummy_sts[sps[0]].lattice,
                           species=["X"]*coords.shape[0],
                           coords = coords,
        )
        st_tmp.merge_sites(mode="del")
        dummy_sts_share.append([st_tmp,sps])

    #for i,st in enumerate(dummy_sts_share):
    #    st.to("POSCAR","{}.dummy.vasp".format(i))

    
    for dir in ALdirs:
        for rpl in range(nreplicas):
            os.chdir(os.path.join(dir, str(rpl)))
            if dir == "AL0":
                energies_ref = np.load("obs_save.npy")[:, 0]
            else:
                energies_ref = np.loadtxt("energy_corr.dat")[:, 1]
            structure_list = [finame for finame in os.listdir() if "structure." in finame]
            step_ids = [int(st_fi.split(".")[1]) for st_fi in structure_list]
            step_ids.sort()
            for i, energy in enumerate(energies_ref):
                structure = Structure.from_file("structure.{}.vasp".format(step_ids[i]))
                mapped_sts = []
                mapping_success = True
                for dummy_st, specs in dummy_sts_share:
                    # perform sublattice by sublattice mapping
                    sp_rm = filter(lambda s: s not in specs, species)
                    st_tmp = structure.copy()
                    st_tmp.remove_species(sp_rm)
                    num_sp = len(st_tmp)
                    # map to perfect lattice for this species
                    st_tmp = map2perflat(dummy_st, st_tmp)
                    st_tmp.remove_species(["X"])
                    mapped_sts.append(st_tmp)
                    if num_sp != len(st_tmp):
                        print("mapping failed for structure {} in replica {}".format(step_ids[i], rpl))
                        mapping_success = False                            
                            
                for sts in mapped_sts[1:]:
                    for i in range(len(sts)):
                        mapped_sts[0].append(
                            sts[i].species_string,
                            sts[i].frac_coords
                        )
                if ignore_species:
                    mapped_sts[0].remove_species(ignore_species)
                if mapping_success:
                    structures.append(mapped_sts[0])
                    energies.append(energy)
            os.chdir(rootdir)

    if len(trainer_input_dirs) == 1:
        generate_inputdir = os.path.join(trainer_input_dirs[0], "generate")
        train_inputdir = os.path.join(trainer_input_dirs[0], "train")
        predict_inputdir = os.path.join(trainer_input_dirs[0], "predict")
    else:
        generate_inputdir = trainer_input_dirs[0]
        train_inputdir = trainer_input_dirs[1]
        predict_inputdir = trainer_input_dirs[2]

    generate_exe = trainer_commands[0]
    train_exe = trainer_commands[1]
    trainer = aenet_trainer(
        structures,
        energies,
        generate_inputdir,
        train_inputdir,
        predict_inputdir,
        generate_exe,
        train_exe,
    )
    trainer.prepare()
    trainer.train()
    trainer.new_baseinput(base_input_dir[0])
    with open("ALloop.progress", "a") as fi:
        fi.write("train\n")
        fi.flush()
        os.fsync(fi.fileno())


def main():
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile)


if __name__ == "__main__":
    main()
