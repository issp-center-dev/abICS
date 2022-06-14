.. highlight:: none

st2abics
-------------------------------

It is sometimes quite tedious to prepare the :ref:`[config] section<config-section>` in abICS :ref:`input files<input_format>`.
To facilitate this, we provide the ``st2abics`` tool, which takes a structure file readable by pymatgen
and converts it to an abICS input template with the ``[config]`` section filled in. An additional control file is required to
tell ``st2abics`` how to break down the original structure file into ``config.base_structure`` and ``config.defect_structure``
(see :ref:`config-section` for definitions). The tool is used as follows::

    $ st2abics -h
    usage: st2abics [-h] inputfi structurefi [outfi]

    Prepare abICS config from structure file

    positional arguments:
      inputfi      toml input file for st2abics
      structurefi  Structure file that can be read by pymatgen Structure.from_file() method
      outfi        Output file to be used as abics input. Defaults to standard output

    optional arguments:
      -h, --help   show this help message and exit

Examples are provided in ``examples/standard/st2abics`` and can be run as follows::

    $ cd examples/standard/st2abics
    $ st2abics st2abics_MgAl2O4.toml MgAl2O4.vasp abics_MgAl2O4.toml # spinel
    $ st2abics st2abics_CuZn.toml CuZn.vasp abics_CuZn.toml # brass
    $ st2abics st2abics_BZY.toml BaZrO3.vasp abics_BZY.toml # Y-doped BaZrO3

The resulting files (abics_MgAl2O4.toml, abics_CuZn.toml, and abics_BZY.toml in the above example) can be used as abICS input after
filling in the ``[mlref]]``, ``[train]``, ``[sampling]`` and ``[observer]`` sections.

Input Format
^^^^^^^^^^^^
Examples of st2abics input files can be found in examples/standard/st2abics
(``st2abics_CuZn.toml``, ``st2abics_MgAl2O4.toml``, and ``st2abics_BZY.toml`` in the above example). 

The format is similar to ``[config]`` section of abICS input file.

Keywords
^^^^^^^^^^
-  ``supercell`` 

   **Format :** list 
   
   **Description :** The size of supercell by the list format [ :math:`\bf{a}, \bf{b}, \bf{c}` ].

-  ``[[config.base_structure]]`` section

   This section specifies how to construct the base_structure that does not exchange
   atoms between lattice sites during the Monte Carlo calculation.

   -  ``species`` 

      **Format :** list of str 
      
      **Description :** Atomic species of the base_structure. The corresponding coordinates
      are extracted automatically from the input structure file. 

   -  ``fix``
   
      **Format :** bool
      
      **Description :** Whether to disallow local relaxation of the base_structure (true) or not (false).

-  ``[[config.defect_structure]]`` section(s)

   This section specifies the sublattice for configurational sampling.
   There can be more than one [[config.defect_structure]] section, e.g., one for cations and one for anions.
  
   -  ``site_center_species``

      **Format :** list of str
      
      **Description :** The species in the original structure file whose coordinates are used as
      lattice sites for configurational sampling.

   -  ``[[config.defect_structure.groups]]`` subsection(s)
      This section specifies the atom groups that reside on the
      lattice sites for configurational sampling. If not provided, it will be constructed automatically from the original
      structure file using ``site_center_species``.

      -  ``name``

         **Format :** str
         
         **Description :** The name of atom group.

      -  ``species``
         
         **Format :** list of str 
         
         **Description :**
         The atom species belonging to the atom group. The default value is a list containing
         only one species specified by ``name``.
         Elements that do not appear in the original structure file can also be specified.
         A vacancy can be represented by an empty list ``[]``.
         As an example, see ``st2abics_BZY.toml`` in the example directory.

      -  ``coords``
      
         **Format :** list of lists of lists or str
         
         **Description :** The coordinates of each atom in the atom group for 
         each orientation that the atom group can take (see description for coords :ref:`here<coords-orr>`). 
         Default value is  ``[[[0.0, 0.0, 0.0]]]``.

      -  ``num``
      
         **Format :** int
         
         **Description :** The number of atom groups of the type specified in this section.
         Make sure to specify the number based on the sites in the supercell, which may be larger than
         the original structure file read in by st2abics.
