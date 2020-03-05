.. highlight:: none

[solver] section
-------------------------------

Specify solver parameters (first-principles calculation) such as solver type (VASP, QE, ...), path to solver, directory with invariant input file.
The example is shown as follows:

  :: 
  
    [solver]
    type = 'vasp'
    path = './vasp'
    base_input_dir = './baseinput'
    perturb = 0.1

Input Format
^^^^^^^^^^^^
Specify a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Keywords
^^^^^^^^^^

    -  ``type``

    **Format :** str

    **Description :**
    Specify the solver type (``OpenMX, QE, VASP``).

    -  ``path``

    **Format :** str

    **Description :**
    Specify the path to the solver.

    -  ``base_input_dir``

    **Format :** str

    **Description :**
    Specify the path to the base input file.

    -  ``perturb``

    **Format :** float

    **Description :**
    If a structure with good symmetry is input, structure optimization tends to stop at the saddle point. In order to avoid this, an initial structure is formed by randomly displacing each atom in proportion to this parameter. It can also be set to 0.0 or false. Default value = 0.0.

