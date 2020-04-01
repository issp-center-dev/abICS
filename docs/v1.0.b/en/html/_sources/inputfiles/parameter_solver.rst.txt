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

       **Format :** str

       **Description :**
       The path to the base input file.

    -  ``perturb``

       **Format :** float

       **Description :**
       If a structure with good symmetry is input, structure optimization tends to stop at the saddle point. In order to avoid this, an initial structure is formed by randomly displacing each atom in proportion to this parameter. It can also be set to 0.0 or false. Default value = 0.0.

