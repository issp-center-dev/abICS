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

import copy
import os
import shlex
import subprocess
import sys
import time
from timeit import default_timer as timer

from mpi4py import MPI

import numpy as np


class runner(object):
    """
    Base class of runner (manager) of external solver program

    Attributes
    ----------
    solver_name : str
        Name of solver program
    path_to_solver : str
        Path to solver program
    run : Any
        Runner object
    nprocs_per_solver : int
        Number of processes which solver program uses
    nthreads_per_proc : int
        Number of threads which a solver process uses
    comm : MPI.Comm
        MPI Communicator
    perturb : float
        perturbation of atom position for structure optimization
    base_solver_input : Solver.Input
        Input manager
    output : Solver.Output
        Output manager
    """

    def __init__(
        self,
        base_input_dir,
        Solver,
        nprocs_per_solver,
        comm,
        perturb=0,
        nthreads_per_proc=1,
        solver_run_scheme="mpi_spawn_ready",
    ):
        """
        Parameters
        -----------
        base_input_dir : str
            Path to the directory including input file templates
        Solver : SolverBase
            Solver
        nprocs_per_solver : int
            Number of processes which one solver program uses
        comm : MPI.Comm
            MPI Communicator
        perturb : float, default 0.0
            Perturbation of atom position
        nthreads_per_proc : int, default 1
            Number of threads which one solver process uses
        solver_run_scheme : str, default "mpi_spawn_ready"
            Scheme how to invoke a solver program

        Raises
        ------
        ValueError
            Raises ValueError if unknown `solver_run_scheme` is passed
        """
        self.solver_name = Solver.name()
        self.path_to_solver = Solver.path_to_solver
        self.base_solver_input = Solver.input
        self.base_solver_input.from_directory(base_input_dir)
        self.nprocs_per_solver = nprocs_per_solver
        self.nthreads_per_proc = nthreads_per_proc
        self.output = Solver.output
        self.comm = comm
        if solver_run_scheme not in Solver.solver_run_schemes():
            print(
                "{scheme} not implemented for {solver}".format(
                    scheme=solver_run_scheme, solver=Solver.name()
                )
            )
            sys.exit(1)
        if solver_run_scheme == "mpi_spawn_ready":
            self.run = run_mpispawn_ready(
                self.path_to_solver, nprocs_per_solver, nthreads_per_proc, comm
            )
        elif solver_run_scheme == "mpi_spawn":
            self.run = run_mpispawn(
                self.path_to_solver, nprocs_per_solver, nthreads_per_proc, comm
            )
        elif solver_run_scheme == "mpi_spawn_wrapper":
            self.run = run_mpispawn_wrapper(
                self.path_to_solver, nprocs_per_solver, nthreads_per_proc, comm
            )
        elif solver_run_scheme == "subprocess":
            self.run = run_subprocess(
                self.path_to_solver, nprocs_per_solver, nthreads_per_proc, comm
            )
        else:
            msg = "Unknown scheme: {}".format(solver_run_scheme)
            raise ValueError(msg)
        self.perturb = perturb

    def submit(self, structure, output_dir):
        """
        Run a solver program and return results

        Parameters
        ----------
        structure : pymatgen.Structure
            Structure of compounds
        output_dir : str
            Name of directory where solver program saves output files

        Returns
        -------
        energy : float
            Total energy
        structure : pymatgen.Structure
            Structure of compounds after optimization
        """
        if self.perturb:
            structure.perturb(self.perturb)
        solverinput = self.base_solver_input
        solverinput.update_info_by_structure(structure)
        self.run.submit(self.solver_name, solverinput, output_dir)
        results = self.output.get_results(output_dir)
        return np.float64(results.energy), results.structure


class runner_multistep(object):
    """
    Sequential runner

    Attributes
    ----------
    runners : list[runner]
        Runners
    """

    def __init__(
        self,
        base_input_dirs,
        Solver,
        runner,
        nprocs_per_solver,
        comm,
        perturb=0,
        nthreads_per_proc=1,
        solver_run_scheme="mpi_spawn_ready",
    ):
        """
        Parameters
        ----------
        base_input_dirs : list[str]
            List of paths to directories including base input files
        Solver : SolverBase
            Solver
        nprocs_per_solver : int
            Number of processes which one solver program uses
        comm : MPI.Comm
            MPI Communicator
        perturb : float, default 0.0
            Perturbation of atom position
        nthreads_per_proc : int, default 1
            Number of threads which one solver process uses
        solver_run_scheme : str, default "mpi_spawn_ready"
            Scheme how to invoke a solver program
        """

        self.runners = []
        assert len(base_input_dirs) > 1
        self.runners.append(
            runner(
                base_input_dirs[0],
                copy.deepcopy(Solver),
                nprocs_per_solver,
                comm,
                perturb,
                nthreads_per_proc,
                solver_run_scheme,
            )
        )
        for i in range(1, len(base_input_dirs)):
            self.runners.append(
                runner(
                    base_input_dirs[i],
                    copy.deepcopy(Solver),
                    nprocs_per_solver,
                    comm,
                    perturb=0,
                    nthreads_per_proc=nthreads_per_proc,
                    solver_run_scheme=solver_run_scheme,
                )
            )

    def submit(self, structure, output_dir):
        energy, newstructure = self.runners[0].submit(structure, output_dir)
        for i in range(1, len(self.runners)):
            energy, newstructure = self.runners[i].submit(
                newstructure, output_dir
            )
        return energy, newstructure

class run_mpispawn:
    """
    Invoker via mpi_comm_spawn

    Attributes
    ----------
    path_to_solver : str
        Path to solver program
    nprocs : int
        Number of process which one solver uses
    nthreads : int
        Number of threads which one solver process uses
    comm : MPI.Comm
        MPI Communicator
    commsize : int
        Size of comm
    commrank : int
        My rank in comm
    worldrank : int
        My rank in MPI.COMM_WORLD
    """

    def __init__(self, path_to_solver, nprocs, nthreads, comm):
        """
        Parameters
        ----------
        path_to_solver : str
            Path to solver program
        nprocs : int
            Number of process which one solver uses
        nthreads : int
            Number of threads which one solver process uses
        comm : MPI.Comm
            MPI Communicator
        """
        self.path_to_solver = path_to_solver
        self.nprocs = nprocs
        self.nthreads = nthreads
        self.comm = comm
        self.commsize = comm.Get_size()
        self.commrank = comm.Get_rank()
        commworld = MPI.COMM_WORLD
        self.worldrank = commworld.Get_rank()

    def submit(self, solver_name, solverinput, output_dir, rerun=2):
        """
        Run solver

        Parameters
        ----------
        solver_name : str
            Name of solver (e.g., VASP)
        solverinput : Solver.Input
            Input manager
        output_dir : str
            Path to directory where a solver saves output
        rerun : int, default = 2
            How many times to restart solver on failed

        Returns
        -------
        status : int
            Always returns 0
        """
        solverinput.write_input(output_dir)

        # Barrier so that spawn is atomic between processes.
        # This is to make sure that vasp processes are spawned one by one according to
        # MPI policy (hopefully on adjacent nodes)
        # (might be MPI implementation dependent...)

        # for i in range(self.commsize):
        #    self.comm.Barrier()
        #    if i == self.commrank:
        failed_dir = []
        cl_argslist = self.comm.gather(
            solverinput.cl_args(self.nprocs, self.nthreads, output_dir), root=0
        )
        solverrundirs = self.comm.gather(output_dir, root=0)

        # checkfilename = "abacus_solver_finished"

        if self.commrank == 0:
            for rundir in solverrundirs:
                solverinput.cleanup(rundir)

            # wrappers = [
            #     "rm -f {checkfile}; {solvername} {cl_args}; echo $? > {checkfile}".format(
            #         checkfile=shlex.quote(
            #             os.path.join(rundir, checkfilename)
            #         ),
            #         solvername=self.path_to_solver,
            #         cl_args=" ".join(map(shlex.quote, cl_args)),
            #     )
            #     for cl_args in cl_argslist
            # ]
            #
            # start = timer()
            # commspawn = [
            #     MPI.COMM_SELF.Spawn(
            #         os.getenv('SHELL'), args=["-c", wrapper], maxprocs=self.nprocs
            #     )
            #     for wrapper in wrappers
            # ]

            start = timer()
            commspawn = [
                MPI.COMM_SELF.Spawn(
                    self.path_to_solver, args=cl_args, maxprocs=self.nprocs
                )
                for cl_args in cl_argslist
            ]
            end = timer()
            print("rank ", self.worldrank, " took ", end - start, " to spawn")
            sys.stdout.flush()
            start = timer()
            for rundir in solverrundirs:
                while True:
                    # if os.path.exists(os.path.join(rundir, checkfilename)):
                    if solverinput.check_finished(rundir):
                        break
                    time.sleep(1)
            end = timer()
            print(
                "rank ",
                self.worldrank,
                " took ",
                end - start,
                " for " + solver_name + "execution",
            )

            # if len(failed_dir) != 0:
            #     print(
            #         solver_name + " failed in directories: \n " + "\n".join(failed_dir)
            #     )
            #     sys.stdout.flush()
            #     if rerun == 0:
            #         MPI.COMM_WORLD.Abort()
        self.comm.Barrier()

        # Rerun if Solver failed
        # failed_dir = self.comm.bcast(failed_dir, root=0)
        # if len(failed_dir) != 0:
        #     solverinput.update_info_from_files(output_dir, rerun)
        #     rerun -= 1
        #     self.submit(solverinput, output_dir, rerun)

        return 0


class run_mpispawn_ready:
    """
    Invoker via mpi_comm_spawn for solvers which is MPI_Comm_spawn-ready

    Attributes
    ----------
    path_to_solver : str
        Path to solver program
    nprocs : int
        Number of process which one solver uses
    nthreads : int
        Number of threads which one solver process uses
    comm : MPI.Comm
        MPI Communicator
    commsize : int
        Size of comm
    commrank : int
        My rank in comm
    worldrank : int
        My rank in MPI.COMM_WORLD
    """

    def __init__(self, path_to_spawn_ready_solver, nprocs, nthreads, comm):
        """
        Parameters
        ----------
        path_to_solver : str
            Path to solver program
        nprocs : int
            Number of process which one solver uses
        nthreads : int
            Number of threads which one solver process uses
        comm : MPI.Comm
            MPI Communicator
        """
        self.path_to_solver = path_to_spawn_ready_solver
        self.nprocs = nprocs
        self.nthreads = nthreads
        self.comm = comm
        self.commsize = comm.Get_size()
        self.commrank = comm.Get_rank()
        commworld = MPI.COMM_WORLD
        self.worldrank = commworld.Get_rank()

    def submit(self, solver_name, solverinput, output_dir, rerun=2):
        """
        Run solver

        Parameters
        ----------
        solver_name : str
            Name of solver (e.g., VASP)
        solverinput : Solver.Input
            Input manager
        output_dir : str
            Path to directory where a solver saves output
        rerun : int, default = 2
            How many times to restart solver on failed

        Returns
        -------
        status : int
            Always returns 0

        Notes
        -----
        If a solver failed (returned nonzero),
        this calls `MPI_Abort` on `MPI_COMM_WORLD` .
        
        """
        solverinput.write_input(output_dir=output_dir)

        # Barrier so that spawn is atomic between processes.
        # This is to make sure that vasp processes are spawned one by one according to
        # MPI policy (hopefully on adjacent nodes)
        # (might be MPI implementation dependent...)

        # for i in range(self.commsize):
        #    self.comm.Barrier()
        #    if i == self.commrank:
        failed_dir = []
        cl_argslist = self.comm.gather(
            solverinput.cl_args(self.nprocs, self.nthreads, output_dir), root=0
        )
        solverrundirs = self.comm.gather(output_dir, root=0)
        if self.commrank == 0:
            start = timer()
            commspawn = [
                MPI.COMM_SELF.Spawn(
                    self.path_to_solver,  # ex. /home/issp/vasp/vasp.5.3.5/bin/vasp",
                    args=cl_args,
                    maxprocs=self.nprocs,
                )
                for cl_args in cl_argslist
            ]
            end = timer()
            print("rank ", self.worldrank, " took ", end - start, " to spawn")
            sys.stdout.flush()
            start = timer()
            exitcode = np.array(0, dtype=np.intc)
            i = 0
            for comm in commspawn:
                comm.Bcast([exitcode, MPI.INT], root=0)
                comm.Disconnect()
                if exitcode != 0:
                    failed_dir.append(solverrundirs[i])
                i = i + 1
            end = timer()
            print(
                "rank ",
                self.worldrank,
                " took ",
                end - start,
                " for " + solver_name + "execution",
            )

            if len(failed_dir) != 0:
                print(
                    solver_name + " failed in directories: \n " + "\n".join(failed_dir)
                )
                sys.stdout.flush()
                if rerun == 0:
                    MPI.COMM_WORLD.Abort()
        self.comm.Barrier()

        # Rerun if Solver failed
        failed_dir = self.comm.bcast(failed_dir, root=0)
        if len(failed_dir) != 0:
            solverinput.update_info_from_files(output_dir, rerun)
            rerun -= 1
            self.submit(solver_name, solverinput, output_dir, rerun)

        # commspawn = MPI.COMM_SELF.Spawn(self.path_to_vasp, #/home/issp/vasp/vasp.5.3.5/bin/vasp",
        #                                args=[output_dir],
        #                                   maxprocs=self.nprocs)

        # Spawn is too slow, can't afford to make it atomic
        # commspawn = MPI.COMM_SELF.Spawn(self.path_to_vasp, #/home/issp/vasp/vasp.5.3.5/bin/vasp",
        #                               args=[output_dir,],
        #                               maxprocs=self.nprocs)
        #        sendbuffer = create_string_buffer(output_dir.encode('utf-8'),255)
        #        commspawn.Bcast([sendbuffer, 255, MPI.CHAR], root=MPI.ROOT)
        # commspawn.Barrier()
        # commspawn.Disconnect()
        # os.chdir(cwd)
        return 0



class run_mpispawn_wrapper:
    """
    Invoker via mpi_comm_spawn

    Attributes
    ----------
    path_to_solver : str
        Path to solver program
    nprocs : int
        Number of process which one solver uses
    nthreads : int
        Number of threads which one solver process uses
    comm : MPI.Comm
        MPI Communicator
    commsize : int
        Size of comm
    commrank : int
        My rank in comm
    worldrank : int
        My rank in MPI.COMM_WORLD
    """

    def __init__(self, path_to_solver, nprocs, nthreads, comm):
        """
        Parameters
        ----------
        path_to_solver : str
            Path to solver program
        nprocs : int
            Number of process which one solver uses
        nthreads : int
            Number of threads which one solver process uses
        comm : MPI.Comm
            MPI Communicator
        """
        self.path_to_solver = path_to_solver
        self.nprocs = nprocs
        self.nthreads = nthreads
        self.comm = comm
        self.commsize = comm.Get_size()
        self.commrank = comm.Get_rank()
        commworld = MPI.COMM_WORLD
        self.worldrank = commworld.Get_rank()

    def submit(self, solver_name, solverinput, output_dir, rerun=2):
        """
        Run solver

        Parameters
        ----------
        solver_name : str
            Name of solver (e.g., VASP)
        solverinput : Solver.Input
            Input manager
        output_dir : str
            Path to directory where a solver saves output
        rerun : int, default = 2
            How many times to restart solver on failed

        Returns
        -------
        status : int
            Always returns 0
        """
        solverinput.write_input(output_dir)

        # Barrier so that spawn is atomic between processes.
        # This is to make sure that vasp processes are spawned one by one according to
        # MPI policy (hopefully on adjacent nodes)
        # (might be MPI implementation dependent...)

        # for i in range(self.commsize):
        #    self.comm.Barrier()
        #    if i == self.commrank:
        failed_dir = []
        cl_argslist = self.comm.gather(
            solverinput.cl_args(self.nprocs, self.nthreads, output_dir), root=0
        )
        solverrundirs = self.comm.gather(output_dir, root=0)

        checkfilename = "abacus_solver_finished"

        if self.commrank == 0:

            wrappers = [
                "rm -f {checkfile}; {solvername} {cl_args}; echo $? > {checkfile}".format(
                    checkfile=shlex.quote(
                        os.path.join(output_dir, checkfilename)
                    ),
                    solvername=self.path_to_solver,
                    cl_args=" ".join(map(shlex.quote, cl_args)),
                )
                for cl_args in cl_argslist
            ]

            start = timer()
            commspawn = [
                MPI.COMM_SELF.Spawn(
                    os.getenv('SHELL'), args=["-c", wrapper], maxprocs=self.nprocs
                )
                for wrapper in wrappers
            ]
            end = timer()
            print("rank ", self.worldrank, " took ", end - start, " to spawn")
            sys.stdout.flush()
            start = timer()
            for rundir in solverrundirs:
                while True:
                    if os.path.exists(os.path.join(rundir, checkfilename)):
                        break
                    time.sleep(1)
            end = timer()
            print(
                "rank ",
                self.worldrank,
                " took ",
                end - start,
                " for " + solver_name + "execution",
            )

            # if len(failed_dir) != 0:
            #     print(
            #         solver_name + " failed in directories: \n " + "\n".join(failed_dir)
            #     )
            #     sys.stdout.flush()
            #     if rerun == 0:
            #         MPI.COMM_WORLD.Abort()
        self.comm.Barrier()

        # Rerun if Solver failed
        # failed_dir = self.comm.bcast(failed_dir, root=0)
        # if len(failed_dir) != 0:
        #     solverinput.update_info_from_files(output_dir, rerun)
        #     rerun -= 1
        #     self.submit(solverinput, output_dir, rerun)

        return 0


class run_subprocess:
    """
    Invoker via subprocess

    Attributes
    ----------
    path_to_solver : str
        Path to solver program
    nprocs : int
        Number of process which one solver uses
    nthreads : int
        Number of threads which one solver process uses
    """

    def __init__(self, path_to_solver, nprocs, nthreads, comm):
        """
        Parameters
        ----------
        path_to_solver : str
            Path to solver program
        nprocs : int
            Number of process which one solver uses
        nthreads : int
            Number of threads which one solver process uses
        comm : MPI.Comm or NoneType
            Never used
        """
        self.path_to_solver = path_to_solver
        self.nprocs = nprocs
        self.nthreads = nthreads

    def submit(self, solver_name, solverinput, output_dir, rerun=0):
        """
        Run solver

        Parameters
        ----------
        solver_name : str
            Name of solver (e.g., VASP)
        solverinput : Solver.Input
            Input manager
        output_dir : str
            Path to directory where a solver saves output
        rerun : int, default = 2
            How many times to restart solver on failed

        Returns
        -------
        status : int
            Always returns 0

        Raises
        ------
        RuntimeError
            Raises RuntimeError when solver failed.
        """
        solverinput.write_input(output_dir=output_dir)
        cwd = os.getcwd()
        os.chdir(output_dir)
        args = solverinput.cl_args(self.nprocs, self.nthreads, output_dir)
        command = [self.path_to_solver]
        command.extend(args)
        to_rerun = False
        print(' '.join(command))
        with open(os.path.join(output_dir, "stdout"), "w") as fi:
            try:
                # subprocess.run(command, check=True, shell=True)
                subprocess.run(' '.join(command), check=True, shell=True)
                # subprocess.run(command, stdout=fi, stderr=subprocess.STDOUT, check=True)
            except subprocess.CalledProcessError as e:
                if rerun > 0:
                    to_rerun = True
                else:
                    raise
        if to_rerun:
            self.submit(solver_name, solverinput, output_dir, rerun - 1)

        os.chdir(cwd)
        return 0
