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

The state of the loop is read from ``ALloop.progress`` and the ``MC*`` directories, so rerunning ``abics_al`` after an
interruption continues where it stopped without repeating finished stages, and a finished loop is extended by
raising ``max_iterations`` and rerunning.

The file format is as follows.

  ::

        [al]
        max_iterations = 5
        sampling_command = "srun -n {nreplicas} abics_sampling {input}"
        walltime = 86400
        walltime_margin = 1800
        resubmit_command = "sbatch job.sh"

Key words
^^^^^^^^^^

- ``max_iterations``

  **Format :** int (natural number, default: 1)

  **Description :** Number of Monte Carlo runs to complete (``MC0`` to ``MC(max_iterations-1)``).
  Can also be given on the command line as ``--max-iterations N``.

- ``train_command``

  **Format :** string (default: ``"abics_train {input}"``)

  **Description :** Command of the training stage; ``{input}`` is replaced by the input file name.
  Include a launcher when resources must be requested, e.g. ``"srun -n 1 --gpus 1 abics_train {input}"`` for GPU training.

- ``sampling_command``

  **Format :** string (default: ``"mpiexec -n {nreplicas} abics_sampling {input}"``)

  **Description :** Command of the sampling stage; ``{nreplicas}`` is replaced by ``nreplicas`` of the ``[sampling]`` section.
  On Slurm use e.g. ``"srun -n {nreplicas} abics_sampling {input}"``.

- ``stage_timeout``

  **Format :** float (seconds, default: 0)

  **Description :** Time after which a training or sampling child process is killed. 0 means no limit.

- ``walltime``

  **Format :** float (seconds, default: 0)

  **Description :** Wall time available to this run of ``abics_al``. 0 means unlimited.
  Before each stage, if the time the same stage took last time plus ``walltime_margin`` does not fit into the
  remaining time, the stage is not started and ``abics_al`` exits with status 2 (finished stages are kept).
  The stage durations are recorded in ``al_driver/durations.json`` and carried over to the next run.

- ``walltime_margin``

  **Format :** float (seconds, default: 600)

  **Description :** Reserve kept for the post-processing of a stage and the end of the job.

- ``resubmit_command``

  **Format :** string (default: empty)

  **Description :** Command run when ``abics_al`` stops because of the wall time (e.g. ``"sbatch job.sh"``).
  Calling ``abics_al`` from a job script and putting the submission command of that script here makes the loop
  continue automatically.

- ``generate_train_input``

  **Format :** bool (default: true)

  **Description :** With ``[train] type = "aenet"`` and a missing ``base_input_dir``, generate a template of the
  training input with ``abics_gen_aenet_input`` before the first training.

Output
^^^^^^^^^^

- ``al_summary.dat``: the error statistics over all replicas of every step from ``AL1`` on (number of configurations,
  RMSE, MAE, largest error), collected from ``AL<n>/convergence.dat``; use it to decide whether to continue the loop.
- ``al_driver/``: the output of the training and sampling child processes (``train_<n>.log``, ``sampling_<n>.log``),
  the stage durations (``durations.json``) and, once the target is reached, a description of the final model,
  ``model_metadata.json`` (abICS version, hash of the input file, number of iterations and training configurations,
  reference solver and trainer settings, species, temperature range and supercell, error statistics of every step).

The exit status is 0 when the target number of iterations is reached, 2 when the run stopped because of the wall time,
and 1 when a stage failed.
