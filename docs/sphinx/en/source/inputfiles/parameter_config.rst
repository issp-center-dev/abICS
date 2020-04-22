.. _config_section:
.. highlight:: none

[config] section
-------------------------------

This section specifies alloy coordination, etc.
An example is shown as follows:

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
Keywords and their values are specified by a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Key words
^^^^^^^^^^

- Specify lattice

    -  ``unitcell``

       **Format :** list

       **Description :**
       Lattice vector :math:`\bf{a}, \bf{b}, \bf{c}` by
       the list format [ :math:`\bf{a}, \bf{b}, \bf{c}` ] .

    -  ``supercell``

       **Format :** list

       **Description :**
       The size of super lattice by the list format[ :math:`\bf{a}, \bf{b}, \bf{c}` ].

- [[config.base_structure]] section

  ``type`` and ``coords`` specify the atomic species that do not move in Monte Carlo calculation and their coordinates.
    If there are multiple atomic species, specify multiple [[config.base_structure]] sections.

    - ``type``

      **Format :** str

      **Description :**  Atomic specie.

    - ``coords``

      **Format :** list of lists or str

      **Description :**
      Coordinates. Specify a list of N elements (number of atoms) arranged in 3 elements representing 3D coordinates, or a string of coordinates arranged in N rows and 3 columns.


- [[config.defect_structure]] section

    This sections specifies the lattice coordinates (coords) and atoms (or atom groups) (groups) that can reside on those lattice sites. Monte Carlo sampling is performed on the lattice specified in this section. In Ver. 1.0, conversion tools from POSCAR and cif will be available.
  
    - ``coords``

      **Format :** list of lists or str

      **Description :**  The coordinates of the lattice sites where atoms reside.
      A list of N elements (number of atoms) arranged in 3 elements representing 3D coordinates, or a string of coordinates arranged in N rows and 3 columns.

    - [[config.defect_structure.groups]] section

      The atom group information to be updated by Monte Carlo.

      -  ``name``

         **Format :** str

         **Description :**
         The name of atomic group.


      -  ``species``

         **Format :** list

         **Description :**
	 The atomic species belonging to the atom group. The default value is a list containing only one specified by ``name``.

      -  ``coords``

	 **Format :** list of lists or str

         **Description :**  The coordinates of each atom in the atom group.
         A list of N elements (number of atoms) arranged in 3 elements representing 3D coordinates, or a string of coordinates arranged in N rows and 3 columns.
	 Default value is  ``[[0.0, 0.0, 0.0]]``.

      - ``relaxation``

	**Format :** list of lists or str

	**Description :**  Whether to optimize structure (coordinates) or not for each atom and dimension.
        A list of N elements (number of atoms) with 3 booleans ("true" or "false"), or a string of "true" or "false" arranged in N rows and 3 columns.
        Default is ``["true", "true", "true"]`` for all the atoms.

      - ``magnetization``

	**Format :** list
	
	**Description :**  Magnetization (the difference between the number of up and down electrons) for each atom.
        Default is 0.0 for all the atoms.
  

      -  ``num``

         **Format :** int

         **Description :**
         The number of atom groups of the type specified in this section.
