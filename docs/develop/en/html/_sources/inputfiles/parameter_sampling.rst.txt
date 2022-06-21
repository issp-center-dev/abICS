.. highlight:: none

[sampling] section
-------------------------------

Specify the parameters of the replica exchange Monte Carlo part, such as the number of replicas, the temperature range, and the number of Monte Carlo steps.
The example is shown as follows.

  ::
  
        [sampling]
        nreplicas = 3
        nprocs_per_replica = 1
        kTstart = 500.0
        kTend = 1500.0
        nsteps = 5
        RXtrial_frequency = 2
        sample_frequency = 1
        print_frequency = 1

Input Format
^^^^^^^^^^^^
Specify a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Keywords
^^^^^^^^^^

- About temperatures

   - ``kTstart``

       **Format :** float (>0)

       **Description :**
       Minimum temperature for the replica.

   - ``kTend``

       **Format :** float (>0)

       **Description :**
       Maximum temperature for the replica.

- About replica 

    - ``nprocs_per_replica``

       **Format :** int (natural number)

       **Description :** The number of processes for the replica. Default value = 1.

    - ``nreplicas``

       **Format :** int (natural number)

       **Description :** The number of replicas.


- Others

   - ``nsteps``

       **Format :** int (natural number)

       **Description :** Number of Monte Carlo steps.

  
   - ``RXtrial_frequency``

       **Format :** int (natural number)

       **Description :** The interval for performing replica exchange trials. For example, setting this value to 1 means that replica exchange is attempted at every Monte Carlo step, while setting this to 2 means that exchange is attempted at every second step. Default = 1.


   - ``sample_frequency``

       **Format :** int (natural number)

       **Description :**     The interval for observation of physical quantities. Default value = 1.

   - ``print_frequency``

       **Format :** int (natural number)

       **Description :**     The interval for saving physical quantities. Default value = 1.

   - ``reload``

       **Format :** bool ("true" or "false")

       **Description :**     Whether to restart a prior calculation from the last step finished last time. Default value = false.


[sampling.solver] section
-------------------------------

This section specifies solver parameters such as solver type (VASP, QE, ...), path to solver, directory with solver-specific input file(s).
An example is shown as follows:

  :: 
  
    [sampling.solver]
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
       The solver type (``OpenMX, QE, VASP, aenet``).

    -  ``path``

       **Format :** str

       **Description :**
       The path to the solver.

    -  ``base_input_dir``

       **Format :** str or list of str

       **Description :**
       The path to the base input file.
       If multiple calculations are set up in the form of a list, each calculation using each input is performed in turn. For the second and subsequent calculations, the structure from the last step of the previous calculation is used as the initial coordinates, and the energy from the last calculation is used. For example, it is possible to perform a fast structural optimization in the first input file at the expense of accuracy, and then perform the structural optimization in the second and later input files with a higher accuracy setting. Or, in the case of grid vector relaxation, one can run the same input multiple times to reset the computational mesh based on a set plane-wave cutoff.

    -  ``perturb``

       **Format :** float

       **Description :**
       If a structure with good symmetry is input, structure optimization tends to stop at the saddle point. In order to avoid this, an initial structure is formed by randomly displacing each atom in proportion to this parameter. It can also be set to 0.0 or false. Default value = 0.0.


    - ``ignore_species``

       **Format :** list

       **Description :**
       Specify atomic species to "ignore" in neural network models such as ``aenet``. For those that always have an occupancy of 1, it is computationally more efficient to ignore their presence when training and evaluating neural network models.

      
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
       Only the specified elements will be passed to ``pw.x`` as command-line options.


