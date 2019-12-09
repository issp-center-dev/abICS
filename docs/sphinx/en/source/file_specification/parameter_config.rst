.. highlight:: none

[config] section
-------------------------------

Specify alloy coordination, etc.
The example is shown as follows:

  ::

    [config]
    unitcell = [[8.1135997772, 0.0000000000, 0.0000000000],
                [0.0000000000, 8.1135997772, 0.0000000000],
                [0.0000000000, 0.0000000000, 8.1135997772]]
    supercell = [1,1,1]

    [[config.base_structure]]
    type = "O"
    coords = [
        [0.237399980, 0.237399980, 0.237399980],
        [0.762599945, 0.762599945, 0.762599945],
        ...
        [0.262599975, 0.262599975, 0.762599945],
        ]

    [[config.defect_structure]]
    coords = [
        [0.000000000, 0.000000000, 0.000000000],
        [0.749999940, 0.249999985, 0.499999970],
        ...
        [0.124999993, 0.624999940, 0.124999993],
        ]
    [[config.defect_structure.groups]]
    name = 'Al'
    # species = ['Al']    # default
    # coords = [[[0,0,0]]]  # default
    num = 16
    [[config.defect_structure.groups]]
    name = 'Mg'
    # species = ['Mg']    # default
    # coords = [[[0,0,0]]]  # default
    num = 8

Input Format
^^^^^^^^^^^^
Specify a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Key words
^^^^^^^^^^

- Specify lattice

    -  ``unitcell``

    **Format :** list

    **Description :**
    Specify lattice vector :math:`\bf{a}, \bf{b}, \bf{c}` by
    the list format [ :math:`\bf{a}, \bf{b}, \bf{c}` ] .

    -  ``supercell``

    **Format :** list

    **Description :**
    Specify the size of super lattice by the list format[ :math:`\bf{a}, \bf{b}, \bf{c}` ].

- [[config.base_structure]] section

  ``type`` and ``coords`` specify the atomic species that do not move in Monte Carlo calculation and their coordinates.
    If there are multiple atomic species, specify multiple [[config.base_strucure]] sections.

    - ``type``

    **Format :** str

    **Description :**  Specify atomic specie.

    - ``coords``

    **Format :** list of lists or str

    **Description :**
    Specify coordinates. Specify a list of N elements (number of atoms) arranged in 3 elements representing 3D coordinates, or a string of coordinates arranged in N rows and 3 columns.


- [[config.defect_structure]] section

    Specify the coordinates (coords) and atoms (group) that can enter the atoms to be updated in Monte Carlo.
    In Ver. 1.0, conversion tools from POSCAR and cif will be available.
  
    - ``coords``

    **Format :** list of lists or str

    **Description :**  Specify the coordinates where the atom enter.
    Specify a list of N elements (number of atoms) arranged in 3 elements representing 3D coordinates, or a string of coordinates arranged in N rows and 3 columns.

    - [[config.defect_structure.groups]] section

        Specify the atom group information to be updated by Monte Carlo.

        -  ``name``

            **Format :** str

            **Description :**
            Specify the name of atomic group.


        -  ``species``

            **Format :** list

            **Description :**
	    Specify the atomic species belonging to the atom group. The default value is a list containing only one specified by ``name``.

        -  ``coords``

            **Format :** list of lists or str

            **Description :**  Specify the coordinates of each atom in the atom group.
            Specify a list of N elements (number of atoms) arranged in 3 elements representing 3D coordinates, or a string of coordinates arranged in N rows and 3 columns.
	    Default value is  `[[0.0, 0.0, 0.0]]`.

        -  ``num``

            **Format :** int

            **Description :**
            Specify the number of this atom group.
