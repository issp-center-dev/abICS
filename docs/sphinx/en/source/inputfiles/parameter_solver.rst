.. highlight:: none

[solver] section
-------------------------------

This section specifies solver parameters such as solver type (VASP, QE, ...), path to solver, directory with solver-specific input file(s).
An example is shown as follows:

  :: 
  
    [solver]
    type = 'vasp'
    path = './vasp'
    base_input_dir = './baseinput'
    perturb = 0.1
    run_scheme = 'mpi_spawn_ready'

Input Format
^^^^^^^^^^^^
Keywords and their values are specified by a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Keywords
^^^^^^^^^^

    -  ``type``

       **Format :** str

       **Description :**
       The solver type (``OpenMX, QE, VASP``).

    -  ``path``

       **Format :** str

       **Description :**
       The path to the solver.

    -  ``base_input_dir``

       **Format :** str or list of str

       **Description :**
       The path to the base input file.
       When more than one directories are specified, they will be executed in sequence during one Monte Carlo step.

    -  ``perturb``

       **Format :** float

       **Description :**
       If a structure with good symmetry is input, structure optimization tends to stop at the saddle point. In order to avoid this, an initial structure is formed by randomly displacing each atom in proportion to this parameter. It can also be set to 0.0 or false. Default value = 0.0.

    - ``run_scheme``

       **Format :** str

       **Description :**
       Way to invoke the solver program.
       For details, please see :ref:`solver_specific_notes`

    -  ``parallel_level`` (Only for QuantumESPRESSO)

       **Format :** dict

       **Description :** 
       How to split parallel cpu resources, i.e., `Parallelization levels <https://www.quantum-espresso.org/Doc/user_guide/node18.html>`_ .
       Key names are long-form command-line options (without the leading hyphen), that is, ``nimage``, ``npools``, ``nband``, ``ntg``, and ``ndiag``.
       Values are the number of parallelization.
       Only the specified elements will be passed to `` pw.x`` as command-line options.

