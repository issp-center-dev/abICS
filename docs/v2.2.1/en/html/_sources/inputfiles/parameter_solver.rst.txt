.. highlight:: none

[sampling.solver], [mlref.solver] section
----------------------------------------------------

These sections specify the parameters of solvers, which calculate the energy of a configuration (e.g., atomic positions).
``sampling.solver`` is used while Monte Carlo sampling, and ``mlref.solver`` is used for generating training dataset of a machine learning model.

In the present version, there are two main types of solvers.

- Ab initio solvers

   - Solvers for calculating DFT energy from atomic position.

   - Indeed, abICS uses an external DFT solver package such as VASP.

      - Parameters of the DFT solver such as convergence criteria are specified by using input files of the solver.

      - Solver specific notes (e.g., input filename) are described in :ref:`solver_specific_notes`.

   - Configurations (atomic positions, atomic species, etc.) are specified by the ``[config]`` section.

- Potts model solver

   - A solver for calculating the energy :math:`E = -\sum_{ij} \delta_{\sigma_i, \sigma_j} \quad (\sigma_i = 0, 1, \dots, Q-1)` of a spin configuration on a hyper cubic lattice (:math:`\{\sigma_i\}`, :math:`\sigma_i = 0, 1, \dots, Q-1`).

   - The dimension and lengths of a hyper cubic lattice :math:`L`, and the local degree of freedom of a spin :math:`Q` are specified in the ``[config]`` section.

   - This solver is for the purpose of testing of algorithms and tutorials.

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
   The solver type (``OpenMX, QE, VASP, aenet, aenetPyLammps, nequip, allegro, mlip_3, potts``).
   When ``potts``, the following parameters are not used.

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


-  ``ignore_species``

   **Format :** list

   **Description :**
   Specify atomic species to "ignore" in neural network models such as ``aenet``. For those that always have an occupancy of 1, it is computationally more efficient to ignore their presence when training and evaluating neural network models.

  
-  ``run_scheme`` (Only for ``sampling.solver``)

   **Format :** str

   **Description :**
   Way to invoke the solver program.
   For details, please see :ref:`solver_specific_notes`

-  ``parallel_level`` (Only for ``type = "QE"``)

   **Format :** dict

   **Description :** 
   How to split parallel cpu resources, i.e., `Parallelization levels <https://www.quantum-espresso.org/Doc/user_guide/node18.html>`_ .
   Key names are long-form command-line options (without the leading hyphen), that is, ``nimage``, ``npools``, ``nband``, ``ntg``, and ``ndiag``.
   Values are the number of parallelization.
   Only the specified elements will be passed to ``pw.x`` as command-line options.

-  ``function`` (Only for ``type = "user"``)

   **Format :** str

   **Description :** 
   Specify a user-defined solver function.

   - The solver function must take a ``pymatgen.core.Structure`` object as an argument and return a ``float`` value.
   - When ``.`` is included, it is assumed that the function is defined in a module, and the module will be automatically imported.

      - For example, if ``function = "mypackage.mymodule.myfunction"`` is specified, ``mypackage.mymodule`` is imported and ``myfunction`` is called.
      - The package is searched from the current directory as well as from the directories specified by the ``PYTHONPATH`` environment variable.
