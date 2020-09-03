from abics.mc_mpi import RefParams
from abics.applications.latgas_abinitio_interface.params import DFTParams, TrainerParams
from abics.applications.latgas_abinitio_interface import aenet_trainer
from abics.applications.latgas_abinitio_interface import map2perflat
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)

from pymatgen import Structure
import numpy as np
import os, sys


def main_impl(tomlfile):

    # We need info from MC steps as well as parameters
    # specific to training
    rxparams = RefParams.from_toml(tomlfile)
    nreplicas = rxparams.nreplicas

    DFTparams = DFTParams.from_toml(tomlfile)
    base_input_dir = DFTparams.base_input_dir

    trainerparams = TrainerParams.from_toml(tomlfile)
    ignore_species = trainerparams.ignore_species
    trainer_commands = trainerparams.exe_command
    trainer_type = trainerparams.solver
    trainer_input_dirs = trainerparams.base_input_dir
    vac_map = trainerparams.vac_map

    configparams = DFTConfigParams.from_toml(tomlfile)
    config = defect_config(configparams)
    perf_st = config.structure

    if trainer_type != "aenet":
        print("Unknown trainer: ", trainer_type)
        sys.exit(1)

    ls = os.listdir()
    rootdir = os.getcwd()
    ALdirs = [dir for dir in ls if "AL" in dir]
    structures = []
    energies = []

    # val_map is a list of list [[sp0, vac0], [sp1, vac1], ...]
    if vac_map:
        vac_map = {specie: vacancy for specie, vacancy in vac_map}
    else:
        vac_map = {}

    for dir in ALdirs:
        for rpl in range(nreplicas):
            os.chdir(os.path.join(dir, str(rpl)))
            energies_ref = np.loadtxt("energy_corr.dat")[:, 1]
            for i, energy in enumerate(energies_ref):
                structure = Structure.from_file("structure_rel.{}.vasp".format(i))
                # map to perfect lattice
                structure = map2perflat(perf_st, structure, vac_map)
                if ignore_species:
                    structure.remove_species(ignore_species)
                structures.append(structure)
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


if __name__ == "__main__":
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile)
