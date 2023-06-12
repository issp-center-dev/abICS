.. highlight:: none

[observer] section
-------------------------------

This section specifies the physical quantity to be calculated.
An example is shown as follows:

  :: 

    [observer]
    [observer.similarity]
    reference_structure = './MgAl2O4.vasp'
    ignored_species = ["O"]

    [[observer.solver]]
    name = "magnetization"
    type = 'aenet'
    path= '~/opt/aenet/bin/predict.x_serial'
    base_input_dir = './baseinput_mag'
    perturb = 0.0
    run_scheme = 'subprocess'
    ignore_species = ["O"]

Input Format
^^^^^^^^^^^^^

Keywords and their values are specified by a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Key words
^^^^^^^^^^

- ``[[observer.solver]]``

  This section specifies a physical quantity to be calculated.
  This section can be specified multiple times.
  This section is the same as the ``sampling.solver`` section except for the ``name`` keyword.

  The quantity with the name ``energy`` is automatically calculated by using ``sampling.solver`` .

  - ``name``

    **Format :** str

    **Description :**
    The name of the physical quantity to be calculated.
    After the calculation is completed, the expected value is output as a file named ``<name>.dat``.

- ``[observer.similarity]``

  "Similarity" is a physical quantity that the ratio of the number of atoms of each element in the same place as the reference state.
  After the calculation is completed, the expected value is output as a file named ``similarity_X.dat`` (``X`` is the element symbol).
  When the ``reference_structure`` keyword is not specified, the similarity is not calculated.

  - ``reference_structure``

    **Format :** str

    **Description :**
    Filename of the structure file of the reference state.

  - ``ignored_species``

    **Format :** list

    **Description :**
    The atom species to be ignored when calculating the similarity.
    For example, if you want to ignore the similarity of oxygen atoms, specify ``["O"]``.
