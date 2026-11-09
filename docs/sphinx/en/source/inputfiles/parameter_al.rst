.. highlight:: none

[al] section
-------------------------------

Settings of ``abics_al``, which runs the whole active-learning loop as one program.
``abics_al`` repeats ``abics_mlref`` in the integrated mode (first-principles calculations),
``abics_train`` (training) and ``abics_sampling`` (Monte Carlo) in the order
``AL0 -> train -> MC0 -> AL1 -> ...`` until ``max_iterations`` Monte Carlo runs exist.
The first-principles calculations are run by this process itself through libEnsemble; training and sampling are
started as child processes with the commands given by ``train_command`` and ``sampling_command``.
Like the integrated mode of ``abics_mlref``, ``abics_al`` must be run **as a single process, without mpiexec / srun**.
libEnsemble is required (``pip install libensemble``).

The ``[mlref]`` section must be set up for the integrated mode (``path`` and ``run_scheme = "subprocess"`` in
``[mlref.solver]``, and ``[mlref.integrated]``); ``abics_al`` overrides ``exec_mode`` with ``"integrated"``.

The state of the loop is read from ``ALloop.progress`` and the ``MC*`` directories, so rerunning ``abics_al`` after an
interruption continues where it stopped without repeating finished stages. A sampling run that was interrupted
(``MC<n>`` exists but is not recorded in ``ALloop.progress``) is moved under ``al_driver/`` and repeated.
A finished loop is extended by raising ``max_iterations`` and rerunning. Conversely, with ``max_iterations`` at or
below the number of recorded Monte Carlo runs the driver stops even in the middle of an iteration (after ``AL<n>``
or its training) and does not run the following stages.
``[sampling] reload = true`` cannot be used (an interrupted sampling run is moved aside and repeated).
A lock file ``al_driver/abics_al.lock`` prevents two ``abics_al`` from running in the same directory.
On SIGTERM (e.g. the scheduler's time limit) ``abics_al`` stops the running stage and releases the lock.
If a lock is left behind (e.g. after SIGKILL), a rerun on the same host checks whether the process is alive
and removes a stale lock; on another host this cannot be checked, so make sure that job has ended and remove
the file by hand. An unrecorded ``MC<n>`` that was modified within the last 120 seconds is not moved aside
but reported as possibly still running (e.g. an ``srun`` child that outlived the driver).

The file format is as follows.

  ::

        [al]
        max_iterations = 5
        sampling_command = "srun -n {nprocs} abics_sampling {input}"
        walltime = 86400
        walltime_margin = 1800
        resubmit_command = "sbatch job.sh"

        [al.stage_estimates]
        mlref = 7200
        train = 1800
        sampling = 3600

Key words
^^^^^^^^^^

- ``max_iterations``

  **Format :** int (natural number, default: 1)

  **Description :** Number of Monte Carlo runs to complete (``MC0`` to ``MC(max_iterations-1)``).
  Can also be given on the command line as ``--max-iterations N``.

- ``train_command``

  **Format :** string (default: ``"abics_train {input}"``)

  **Description :** Command of the training stage; ``{input}`` is replaced by the input file name
  (``{nreplicas}`` and ``{nprocs}`` are available too).
  Include a launcher when resources must be requested, e.g. ``"srun -n 1 --gpus 1 abics_train {input}"`` for GPU training
  (if ``[train.exe_command]`` also uses ``srun``, the nested Slurm job step may need ``--overlap``).

- ``sampling_command``

  **Format :** string (default: ``"mpiexec -n {nprocs} abics_sampling {input}"``)

  **Description :** Command of the sampling stage. ``{nprocs}`` is replaced by the number of MPI processes
  ``abics_sampling`` needs (``nreplicas`` of the ``[sampling]`` section times the number of ensemble members,
  normally ``nreplicas``) and ``{nreplicas}`` by ``nreplicas``. On Slurm use e.g.
  ``"srun -n {nprocs} abics_sampling {input}"``.

  Commands are split on whitespace and run without a shell. Only the placeholders above are substituted;
  other braces (e.g. ``--cpu-bind=map_cpu:{0,1}``) are left as they are. ``~`` and environment variables are
  expanded in every argument.

- ``stage_timeout``

  **Format :** float (seconds, default: 0)

  **Description :** Time after which a training or sampling child process is stopped. 0 means no limit.
  The process receives SIGTERM (so that ``mpirun`` / ``srun`` can clean up their ranks) and SIGKILL 30 seconds later
  if it is still alive; the stage then counts as failed (exit status 1).

- ``walltime``

  **Format :** float (seconds, default: 0)

  **Description :** Wall time available to this run of ``abics_al``. 0 means unlimited.
  Before each stage, if the expected duration of the stage plus ``walltime_margin`` does not fit into the
  remaining time, the stage is not started and ``abics_al`` exits with status 2 (finished stages are kept).
  The expected duration is the longest the stage has taken in any earlier run (recorded in
  ``al_driver/durations.json`` and carried over to the next run), or ``stage_estimates`` for a stage that
  has not run yet. Without either, the stage is started with a warning.
  If a stage cannot be started before any stage was finished in this run, resubmitting would never make
  progress, so ``abics_al`` exits with status 1 instead of resubmitting: raise ``walltime``, lower the
  estimate, or forget outdated records with ``abics_al --forget-durations`` (or by removing
  ``al_driver/durations.json``), also after reducing ``nsteps`` or the like.

- ``walltime_margin``

  **Format :** float (seconds, default: 600)

  **Description :** Reserve kept for the post-processing of a stage and the end of the job.

- ``stage_estimates``

  **Format :** table (seconds keyed by ``mlref``, ``train``, ``sampling``; default: empty)

  **Description :** Initial estimates of the stage durations, used for the wall-time decision in the first job
  (e.g. ``sampling = 3600`` in ``[al.stage_estimates]``).

- ``progress_wait``

  **Format :** float (seconds, default: 30)

  **Description :** Maximum time to wait, after a training or sampling child process exited, for the
  ``ALloop.progress`` it wrote to become visible (delay of shared file systems).

- ``resubmit_command``

  **Format :** string (default: empty)

  **Description :** Command run when ``abics_al`` stops because of the wall time (e.g. ``"sbatch job.sh"``).
  Calling ``abics_al`` from a job script and putting the submission command of that script here makes the loop
  continue automatically. If this command fails (non-zero status), ``abics_al`` exits with status 1.

- ``generate_train_input``

  **Format :** bool (default: true)

  **Description :** With ``[train] type = "aenet"`` and a missing (or empty) ``base_input_dir``, generate a template
  of the training input with ``abics_gen_aenet_input`` before the first training (with a warning). Set it to false
  when you provide the input yourself, so that a misspelt ``base_input_dir`` is not hidden.

Output
^^^^^^^^^^

- ``al_summary.dat``: the error statistics over all replicas of every step from ``AL1`` on (number of configurations,
  RMSE, MAE, largest error), collected from ``AL<n>/convergence.dat``; use it to decide whether to continue the loop.
- ``al_driver/``: the output of the training and sampling child processes (``train_<n>.log``, ``sampling_<n>.log``,
  appended on reruns), the stage durations (``durations.json``), interrupted sampling runs that were moved aside
  (``MC<n>.interrupted-<timestamp>``), the driver's own log (``abics_al_<timestamp>.log``; the driver does not
  write to ``[log] logfile_path`` because the training and sampling programs overwrite it) and, once the target is
  reached, a description of the current model, ``model_metadata.json`` (abICS version, hash of the input file,
  the number of completed iterations ``iterations`` and of AL steps the model was trained on ``trained_on_steps``
  (larger than ``iterations`` when the loop was stopped in the middle of an iteration), the number of reference
  calculations, reference solver and trainer settings, species, temperature range and supercell, error statistics
  of every step). Note that the statistics of ``AL<n>`` evaluate the model trained after ``AL<n-1>`` and used in
  ``MC<n-1>``; the last trained model itself is not evaluated.

The exit status is 0 when the target number of iterations is reached, 2 when the run stopped because of the wall time,
and 1 when a stage failed (or the input is wrong, or the resubmission command failed).
