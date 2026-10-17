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

"""
``abics_mlref``: the reference-data step of the staged active-learning workflow.

Each invocation performs one stage and exits:

* no ``MC*`` and no ``AL0``: generate random configurations and write the
  solver inputs of ``AL0``;
* ``AL0`` exists but ``MC0`` does not: parse the solver outputs of ``AL0``
  (or write the inputs of the next relaxation stage);
* ``MC<n-1>`` is the newest MC run and ``AL<n>`` does not exist: select
  configurations from ``MC<n-1>`` and write the solver inputs of ``AL<n>``;
* ``AL<n>`` exists: parse its solver outputs (or write the next stage).

The current stage is inferred from the directories and progress files in
the working directory.  The per-configuration work is done by
:mod:`abics.applications.latgas_abinitio_interface.activelearn`.
"""

from typing import MutableMapping

import sys
import os
import datetime

import numpy as np

from abics import __version__

from abics.mc_mpi import (
    RX_MPI_init,
    RefParams,
)
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)
from abics.applications.latgas_abinitio_interface.base_solver import SolverBase, create_solver
from abics.applications.latgas_abinitio_interface.params import ALParams, DFTParams, TrainerParams
from abics.applications.latgas_abinitio_interface import activelearn as al

from abics.util import exists_on_all_nodes
from pymatgen.core import Structure

import logging
import abics.loggers as loggers

logger = logging.getLogger("main")


def newest_mc_index(comm) -> int:
    """Index of the MC run that does not exist yet (``MC0``, ``MC1``, ...)."""
    index = 0
    while exists_on_all_nodes(comm, f"MC{index}"):
        index += 1
    return index


def finish(comm, replica: int, rootdir: str, progress_label: str, messages) -> None:
    """Record ``progress_label`` in ``ALloop.progress`` on rank 0 and exit."""
    comm.Barrier()
    if replica == 0:
        logger.info("-Writing ALloop.progress")
        al.record_progress(rootdir, progress_label)
        for msg in messages:
            logger.info(msg)
        now = datetime.datetime.now()
        logger.info(f"Exiting normally on {now}.\n")
    sys.exit(0)


def finish_stage(comm, replica: int, aldir: str, stage: int, messages) -> None:
    """Record the finished relaxation stage in ``baseinput.progress`` on rank 0 and exit."""
    comm.Barrier()
    if replica == 0:
        al.record_stage(aldir, stage)
        for msg in messages:
            logger.info(msg)
        now = datetime.datetime.now()
        logger.info(f"Exiting normally on {now}.\n")
    sys.exit(0)


def finish_inputs(comm, replica: int, rootdir: str, rundir_list, messages) -> None:
    """Gather the solver input directories to ``rundirs.txt`` on rank 0 and exit."""
    rundir_list = comm.gather(rundir_list, root=0)
    if replica == 0:
        rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
        al.write_rundirs(rootdir, rundir_list)
        for msg in messages:
            logger.info(msg)
        now = datetime.datetime.now()
        logger.info(f"Exiting normally on {now}.\n")
    sys.exit(0)


def main_impl(params_root: MutableMapping):
    rxparams = RefParams.from_dict(params_root["mlref"])
    ndata = rxparams.ndata
    # the sampler that produced the MC* directories (same default as abics_sampling)
    sampler_type = params_root["sampling"].get("sampler", "RXMC")
    comm = RX_MPI_init(rxparams.nreplicas, rxparams.seed)
    alparams = ALParams.from_dict(params_root["mlref"]["solver"])
    mcparams = DFTParams.from_dict(params_root["sampling"]["solver"])
    configparams = DFTConfigParams.from_dict(params_root["config"])
    trainerparams = TrainerParams.from_dict(params_root["train"])

    energy_ref_empty = trainerparams.energy_ref
    mult_now = np.prod(configparams.supercell)

    solver: SolverBase = create_solver(alparams.solver, alparams)

    perturb = alparams.perturb
    solver_output = solver.output
    solver_input = solver.input
    nstages = len(alparams.base_input_dir)

    myreplica = comm.Get_rank()

    # Find newest MC run
    nextMC_index = newest_mc_index(comm)

    rootdir = os.getcwd()
    if nextMC_index == 0:  # "MC0" does not exist
        # Random sampling!
        fmtstr = al.input_dir_format(ndata)
        aldir = os.path.join(rootdir, "AL0")
        replica_dir = os.path.join(aldir, str(myreplica))

        if exists_on_all_nodes(comm, "AL0"):
            if exists_on_all_nodes(comm, "ALloop.progress"):
                logger.error("It seems you've already run the first active learning step. You should train now.")
                sys.exit(1)

            # attempting post processing of ALed output
            # (wait until baseinput.progress is seen consistently from all ranks)
            exists_on_all_nodes(comm, os.path.join(aldir, al.STAGE_FILE))
            runstep = al.read_stage(aldir)
            finalrun = runstep + 1 == nstages  # we're done with this AL step
            if not finalrun:
                solver_input.from_directory(alparams.base_input_dir[runstep + 1])

            logger.info(f"-Parsing {alparams.solver} results in AL0/*/input*/baseinput{runstep}...")
            if finalrun:
                logger.info(f"--This is the final {alparams.solver} calculation step for AL0")
            else:
                logger.info("--Input files for the next calculation step will be")
                logger.info(f"---prepared in AL{nextMC_index}/*/input*/baseinput{runstep+1}")

            energies = []
            rundir_list = []
            for i in range(ndata):
                input_dir = os.path.join(replica_dir, fmtstr.format(i))
                if finalrun:
                    # energy_calculator may return the structure w/o ignored species
                    energy, st = al.parse_reference_output(
                        solver_output, input_dir, runstep, alparams.ignore_species
                    )
                    st.sort(key=lambda site: site.species_string)
                    st.to(
                        fmt="POSCAR",
                        filename=os.path.join(replica_dir, f"structure.{i}.vasp"),
                    )
                    energies.append(energy)
                else:
                    energy, st = al.parse_reference_output(solver_output, input_dir, runstep)
                    rundir_list.append(
                        al.write_next_stage_input(solver_input, st, input_dir, runstep + 1)
                    )

            if finalrun:
                al.write_energy_corr(
                    os.path.join(replica_dir, "energy_corr.dat"),
                    [("NaN", energy, i) for i, energy in enumerate(energies)],
                )
                finish(
                    comm, myreplica, rootdir, "AL0",
                    ["--Done. Please run abics_train next."],
                )
            else:
                rundir_list = comm.gather(rundir_list, root=0)
                if myreplica == 0:
                    rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                    al.write_rundirs(rootdir, rundir_list)
                finish_stage(
                    comm, myreplica, aldir, runstep,
                    [
                        f"-Finished preparing {alparams.solver} input in AL0/*/input*/baseinput{runstep+1}.",
                        "--See rundirs.txt for a full list of the directories.",
                        "-Please perform the calculations in those directories before running abics_mlref again.",
                    ],
                )
        else:
            # "AL0" nor "MC0" does not exist:
            # this should be the first time running this script
            if myreplica == 0:
                logger.info("-No preceding MC or AL runs exist. Preparing {} input by random sampling".format(alparams.solver))
                logger.info("--Input files will be prepared in AL0/*/input*/baseinput0")

            config = defect_config(configparams)
            comm.Barrier()
            if myreplica == 0:
                os.mkdir(aldir)

            # Wait until the new directory is visible from all ranks
            if not exists_on_all_nodes(comm, "AL0"):
                logger.error("Failed to create AL0 directory.")
                sys.exit(1)

            solver_input.from_directory(alparams.base_input_dir[0])
            rundir_list = []
            for i, structure in enumerate(al.generate_random_structures(config, ndata)):
                input_dir = os.path.join(replica_dir, fmtstr.format(i))
                rundir_list.append(
                    al.write_reference_input(solver_input, structure, perturb, input_dir)
                )
            finish_inputs(
                comm, myreplica, rootdir, rundir_list,
                [
                    f"-Finished preparing {alparams.solver} input in AL0/*/input*/baseinput0.",
                    "--See rundirs.txt for a full list of the directories.",
                    "--Please perform the calculations in those directories before running abics_mlref again.",
                ],
            )

    else:  # "MC*" exists
        # Active learning!
        MCdirname = f"MC{nextMC_index-1}"
        mc_replica_dir = os.path.join(rootdir, MCdirname, str(myreplica))
        last_li = al.last_progress(rootdir)
        if last_li != MCdirname:
            logger.error("The last action:")
            logger.error(f"  expected: {MCdirname}")
            logger.error(f"  actual:   {last_li}")
            logger.error("You shouldn't run activelearn now.")
            logger.error("Either train (abics_train) or MC (abics_sampling) first.")
            sys.exit(1)
        # Configurations sampled by MC{n-1} in this replica, independent of the sampler
        candidates = al.reference_candidates(mc_replica_dir, sampler_type)
        nsamples = len(candidates)
        fmtstr = al.input_dir_format(nsamples)
        # energy predicted by the current model, shifted to the solver's energy reference
        energy_ref = {
            candidate.index: candidate.energy + energy_ref_empty * mult_now
            for candidate in candidates
        }
        ALstep = nextMC_index
        aldir = os.path.join(rootdir, f"AL{ALstep}")
        replica_dir = os.path.join(aldir, str(myreplica))
        # NOTE: the Config object is not used below, but constructing it calls
        # Config.shuffle() and thus advances the random number generator.
        # It is kept so that the random perturbations (and the "random"
        # sampler) stay identical to those of earlier versions of abics_mlref.
        config = defect_config(configparams)

        if exists_on_all_nodes(comm, replica_dir):  # We're going to read solver results
            # (wait until baseinput.progress is seen consistently from all ranks)
            exists_on_all_nodes(comm, os.path.join(aldir, al.STAGE_FILE))
            runstep = al.read_stage(aldir)
            finalrun = runstep + 1 == nstages  # This is the last run for this AL step
            if not finalrun:
                solver_input.from_directory(alparams.base_input_dir[runstep + 1])

            logger.info(f"-Parsing {alparams.solver} results in AL{nextMC_index}/*/input*/baseinput{runstep}...")
            if finalrun:
                logger.info(f"--This is the final {alparams.solver} calculation step for AL{nextMC_index}")
            else:
                logger.info("--Input files for the next calculation step will be")
                logger.info(f"---prepared in AL{nextMC_index}/*/input*/baseinput{runstep+1}")

            energy_corrlist = []
            relax_max = []
            rundir_list = []
            sample_list = al.list_input_indices(replica_dir)
            missing = [i for i in sample_list if i not in energy_ref]
            if missing:
                logger.error(
                    f"No sample of {MCdirname}/{myreplica} corresponds to the input directories "
                    f"{[fmtstr.format(i) for i in missing]} of AL{ALstep}/{myreplica}."
                )
                logger.error("The MC output or the sampler setting ([sampling] sampler) has changed since the inputs were generated.")
                comm.Abort(1)
            for i in sample_list:
                input_dir = os.path.join(replica_dir, fmtstr.format(i))
                if finalrun:
                    # structure before perturbation: used to restore the ignored species
                    # and to measure the relaxation
                    st_in = Structure.from_file(os.path.join(input_dir, "initial.vasp"))
                    # energy_calculator may return the structure w/o ignored species
                    energy, st_rel = al.parse_reference_output(
                        solver_output, input_dir, runstep, alparams.ignore_species, st_in
                    )
                    st_rel.to(
                        fmt="POSCAR",
                        filename=os.path.join(replica_dir, f"structure.{i}.vasp"),
                    )
                    energy_corrlist.append((energy_ref[i], energy, i))
                    relax_max.append(al.max_displacement(st_in, st_rel))
                else:
                    energy, st_rel = al.parse_reference_output(solver_output, input_dir, runstep)
                    rundir_list.append(
                        al.write_next_stage_input(solver_input, st_rel, input_dir, runstep + 1)
                    )
            rundir_list = comm.gather(rundir_list, root=0)
            if myreplica == 0 and not finalrun:
                rundir_list = [rundir for sublist in rundir_list for rundir in sublist]
                al.write_rundirs(rootdir, rundir_list)

            if finalrun:
                al.write_energy_corr(os.path.join(replica_dir, "energy_corr.dat"), energy_corrlist)
                al.write_relax_max(os.path.join(replica_dir, "relax_max.dat"), relax_max)
                # Errors of the model energies against the solver, over all replicas
                # (gathered in memory: rank 0 must not depend on seeing the files
                # just written by the other ranks on a shared file system)
                replica_rows = comm.gather(energy_corrlist, root=0)
                if myreplica == 0:
                    per_replica, overall = al.al_step_error_statistics(replica_rows)
                    lines = al.format_convergence(f"AL{ALstep}", per_replica, overall)
                    al.write_convergence(os.path.join(aldir, "convergence.dat"), lines)
                    logger.info(f"-Errors of the model energies against {alparams.solver} results (AL{ALstep}/convergence.dat):")
                    for line in lines:
                        logger.info(f"  {line}")
                finish(
                    comm, myreplica, rootdir, f"AL{ALstep}",
                    ["--Done. Please run abics_train next."],
                )
            else:
                finish_stage(
                    comm, myreplica, aldir, runstep,
                    [
                        f"-Finished preparing {alparams.solver} input in AL{nextMC_index}/*/input*/baseinput{runstep+1}.",
                        "--See rundirs.txt for a full list of the directories.",
                        "--Please perform the calculations in those directories before running abics_mlref again.",
                    ],
                )

        else:  # Create first input for this AL step
            logger.info(f"-This is the first run for AL{nextMC_index}")
            logger.info(f"--Configurations from MC{nextMC_index-1} will be converted to {alparams.solver} input")
            os.makedirs(replica_dir, exist_ok=False)
            solver_input.from_directory(alparams.base_input_dir[0])
            selected = al.select_candidates(candidates, rxparams)
            rundir_list = []
            # We read structures from MC, but we also need to add ignored species.
            # We utilize structure_norel.*.vasp, which is the structure before solver run
            # in MC with ignored species.
            structures = al.load_candidate_structures(selected, mcparams.ignore_species)
            for candidate, structure in zip(selected, structures):
                input_dir = os.path.join(replica_dir, fmtstr.format(candidate.index))
                rundir_list.append(
                    al.write_reference_input(solver_input, structure, perturb, input_dir)
                )
            finish_inputs(
                comm, myreplica, rootdir, rundir_list,
                [
                    f"-Finished preparing {alparams.solver} input in AL{nextMC_index}/*/input*/baseinput0.",
                    "--See rundirs.txt for a full list of the directories.",
                    "--Please perform the calculations in those directories before running abics_mlref again.",
                ],
            )


def main():

    now = datetime.datetime.now()

    import toml
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    params_root = toml.load(tomlfile)

    loggers.set_log_handles(
        app_name="activelearn",
        level=logging.INFO,
        params=params_root.get("log", {}))

    logger.info(f"Running abics_mlref (abICS v{__version__}) on {now}")
    logger.info(f"-Reading input from: {tomlfile}")

    main_impl(params_root)


if __name__ == "__main__":
    main()
