.. highlight:: none

Download
~~~~~~~~~~~~~~~~~~~~~~

The source codes of abICS can be obtained from `GitHub page <https://github.com/issp-center-dev/abICS>`_ .

``$ git clone https://github.com/issp-center-dev/abICS``

Prerequisites
~~~~~~~~~~~~~~~~~~~~~~

- python3
- numpy
- scipy
- toml (for parsing input files)
- mpi4py (for parallel tempering)
- pymatgen (for parsing vasp I/O)
- qe-tools (for parsing QE I/O)

To use VASP as a solver, a patch must be applied to use MPI_COMM_SPAWN. If you wish to use it, please contact us (the e-mail address is written in :doc:`../contact/index` ).
  
Directory structure
~~~~~~~~~~~~~~~~~~~~~~

The directory structure of abICS is given as follows:

  :: 

     - examples
        - standard
            - spinel
        - expert 
            - ising2D
            - 2D_hardcore
            …
    - make_wheel.sh
    - abics
        - python module

``examples/standard`` contains samples that can be run by simple files.
``examples/expert`` contains examples by using python module.

A set of python modules are located in the ``py_mc`` directory.


      
Install
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Make wheel file by typing following command:

``$ ./make_wheel.sh``

2. Install using the created file as follows:

``$ pip install dist/abics-\*.whl``

If you want to change the install directory, use
``--user`` option or ``--prefix = DIRECTORY`` ( ``DIRECTORY`` is the path to the directory where you want to install) option. In the following, the case for using ``--user`` option is shown:

``$ pip install --user dist/abics-\*.whl``
