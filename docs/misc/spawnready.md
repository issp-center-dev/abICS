Making legacy fortran/c/c++ code spawn-ready
============================================
Note: some basic knowledge of MPI and mpi4py is assumed in this document.

## Background: spawning parallel programs from parallel programs in parallel

abICS relies on multilayered parallelism in order to take advantage of massively parallel supercomputer systems. The program is started with N<sub>rep</sub> MPI processes that handles the sampling for each replica or walker. Each of those processes, in turn, spawns parallel solvers (DFT, classical MD, etc.) that are also written using MPI to perform the energy calculations necessary for thermodynamic sampling. 

At first sight, the idea seems simple enough; each of the N<sub>rep</sub> processes can call the solvers by invoking the mpirun command (or mpiexec or whatever command used for running mpi programs) using the "system()" function that is provided in most programming languages. For example, in python, you may write an MPI program that calls another MPI program like this:
```python
test_system.py:

from mpi4py import MPI
import os
os.system('mpirun -np 2 some_mpi_program')
```
and run it like this:
```
$ mpirun -np 2 python test_system.py
```
However, whether this works or not is highly dependent on the system configuration and MPI implementation. For example, Open MPI will fail with an error message `Open MPI does not support recursive calls of mpirun`. Even if the implementation allows for recursive calls, it is usually quite nontrivial to specify on which node/processor the MPI program should run in a multi-node scenario.

## Using MPI_Comm_spawn

To allow for such use cases, the MPI standard provides the `MPI_Comm_spawn` function. This function (or subroutine if you are using fortran) creates a new group of child processes with its own MPI_COMM_WORLD and establishes an intercommunicator with the parent communicator. Here's an example using python with mpi4py:

```python
test_spawn.py:

from mpi4py import MPI
commspawn = MPI.COMM_SELF.Spawn(
    "some_parallel_code",
    args = ["some_argument"],
    maxprocs = 4
    )
```
which may be run as 
```
$ mpirun -np 10 python test_spawn.py
```
This will first run with 10 MPI processes, each of which will spawn 4 processes of `some_parallel_code`, meaning that there will be 50 processes in total and 11 `MPI_COMM_WORLD`s. If you replace `MPI.COMM_SELF` with `MPI.COMM_WORLD`, then only 4 processes will be spawned leading to 14 processes in total with 2 `MPI_COMM_WORLD`s; this is because `MPI_Comm_spawn` is collective over the parent communicator. Some MPI implementations require input of the total number of processes that can be spawned. See section "Comments on MPI implementation" in abICS manual for details.

So, does `MPI_Comm_spawn` solve all of our problems? The answer is NO. If one wants to spawn legacy MPI code without any modification, there are mainly three issues, all pertaining to the limited control provided by the `MPI_Comm_spawn` function:
1. `MPI_Comm_spawn` has no way to specify the working directory (at least, the MPI standard does not specify how)
2. `MPI_Comm_spawn` does not provide a way to redirect stdout or stderr
3. `MPI_Comm_spawn` does not provide a way to detect  termination of the child processes

Issue 1 is a complete deal-breaker if your legacy code relies on the current working directory for finding input files. In fact, the MPI standard does not even specify where the working directory of the spawned processes should be. A prominent DFT code that falls in this category is VASP, as it looks in the current working directory for all of it's input files (`POSCAR`, `INCAR`, `KPOINTS`, `POTCAR`) whose names are hard-coded in the program. On the other hand, many other codes allow specification of the input files on the command line; in those cases, the spawning program can specify the absolute path of the input files as command line arguments.

Issue 2 means that the stdout and stderr of the parent and all child processes are jumbled up in the stdout and stderr streams. This is problematic if it is necessary to parse the stdout for extracting information such as the total energy used in Metropolis sampling. If it is OK to throw out the stdout, this is not a problem.

Issue 3 comes from the fact that `MPI_Comm_spawn` returns right after the child processes are spawned.  One way to detect termination of the child processes without modification of the legacy code is to monitor file output at preset intervals. This scheme can only be used if a certain file with certain contents is output right before the program terminates.

If these issues can all be resolved by the above-mentioned workarounds, then `MPI_Comm_spawn` can be used without modification of the legacy MPI code. If not, it will be necessary to make the code "spawn-ready" by modifying the source code as follows. abICS provides two schemes for using `MPI_Comm_spawn` to spawn child processes. One of them, which we will call the "MPI_Spawn" scheme, requires no modification on the legacy code side, but the abovementioned workarounds have to be implemented on the abICS side. The second scheme, which we call the "MPI_Spawn_ready" scheme, requires modification of the legacy code but the implementation of the interface on the abICS side is much simpler and fault-tolerant/foolproof.

## Making legacy MPI code spawn-ready

In the following, we will outline the protocol for making legacy code spawn-ready to be used from abICS in the "MPI_Spawn_ready" scheme; the same protocol should be useful for a variety of other `MPI_Comm_spawn` scenarios that want to utilize legacy MPI programs with minimum modification. The code examples will be in Fortran90, although it should be fairly trivial to write the same code in C/C++.

### Taking the working directory as command line argument
If the legacy code relies on the current working directory to find the input files (like VASP mentioned above), then we need to modify the code to take a command line argument that specifies the working directory and to subsequently change into that directory. This can be realized by adding to the beginning of the main program something like the following code:
```fortran
! In variable declaration section
CHARACTER spawndir*255

! At the beginning of the statements section
if(iargc()>0) then
    call getarg(1, spawndir)
    CALL chdir(TRIM(spawndir))
endif
```
Note that this code snippet works without modification only if the original code takes no command line arguments. In other cases, you need to write appropriate logic to take an additional argument corresponding to the working directory.

### Redirecting stdout and stderr
The following code redirects stdout and stderr to file if the program is spawned.
```fortran
! In declaration section
INTEGER :: mpi_comm_parent ! parent communicator (if spawned)
INTEGER :: IERR

! Near beginning of executable statements
CALL MPI_COMM_GET_PARENT(mpi_comm_parent,IERR)
IF (mpi_comm_parent .NE. MPI_COMM_NULL) THEN
    OPEN(UNIT=6,FILE='stdout',STATUS='UNKNOWN')
    OPEN(UNIT=0,FILE='stderr',STATUS='UNKNOWN')
ENDIF
```
Note that the unit number for stdout is usually 6 and stderr is 0. This means that writes to UNIT * will get redirected to 'stdout'. However, this is not specified by the Fortran standard, so there can be exceptions. In those cases, we need to specify the correct unit number for that environment. Also the default behavior for the Intel Compiler is that UNIT 6 and * are treated differently, which actually goes against the Fortran 2003 standard. To get the proper behavior when using the Intel Compiler, it is necessary to set the compiler option "-assume noold_unit_star".

### Communicating termination of the program
In abICS, we require that the child processes returns an exit code through the parent-child intercommunicator using MPI_BCAST. We also require an MPI_COMM_DISCONNECT call to free up the communicator.
```fortran
! In declaration section
INTEGER id_in_group, mpi_comm_parent, excode

! Before MPI_FINALIZE
! At this point, excode should be set to 0 
! if no errors, and it should be > 0 when 
! we want to report that there was a 
! problem in the execution
CALL MPI_COMM_GET_PARENT(mpi_comm_parent,ierror)
IF (mpi_comm_parent .NE. MPI_COMM_NULL) THEN
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, id_in_group, ierror)
    IF (id_in_group==0) THEN
        CALL MPI_BCAST(excode, 1, MPI_INTEGER, MPI_ROOT, &
            mpi_comm_parent, ierror)
    ELSE
        CALL MPI_BCAST(excode, 1, MPI_INTEGER, MPI_PROC_NULL, &
            mpi_comm_parent, ierror)
    ENDIF
    CALL MPI_COMM_DISCONNECT(mpi_comm_parent,ierror)
ENDIF
```

## Remaining issues and idealistic solutions
We have found that `MPI_Comm_spawn` is actually quite slow (probably just as slow as an mpirun call), especially when multiple processes call the function at the same time. Ideally, our use case should not require `MPI_Comm_spawn`; if the DFT code provides a library interface (i.e., subroutine or function) that takes an MPI intracommunicator as an argument to be used as COMM_WORLD, then abICS can completely avoid using `MPI_Comm_spawn`. In this case, abICS will be run with the maximum number of processes available to the calculation job; it can then group processes into intracommunicator groups and make library calls to the DFT library to perform the calculations without spawning new processes. However, making legacy MPI programs into subroutines that can be called from other programs is actually quite nontrivial. Issues that immediately come to mind include removing all `STOP`, `MPI_INIT`, and `MPI_FINALIZE` statements (a `STOP` statement in the subroutine terminates the calling program entirely), handling of module variables, and deallocation of all dynamically allocated arrays. 

Another possibility is to make the DFT program into a server that takes input, performs calculations on the input, then waits for the next input without terminating. The Qbox code (http://qboxcode.org) seems to provide such an interface.

Our conclusion at this point is that it is easiest to use `MPI_Comm_spawn`, although this may change as DFT codes evolve to keep pace with new computing paradigms (e.g., using interpreter languages like `Python` as glue to integrate multiple programs into custom made workflows).
