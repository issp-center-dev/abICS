.. highlight:: none

[replica] section
-------------------------------

Specify the parameters of the replica exchange Monte Carlo part, such as the number of replicas, the temperature range, and the number of Monte Carlo steps.
The example is shown as follows.

  ::
  
        [replica]
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

    -  ``kTstart``

    **Format :** float (>0)

    **Description :**
    Minimum temperature for the replica.

    -  ``kTend``

    **Format :** float (>0)

    **Description :**
    Maximum temperature for the replica.

    -  ``nsteps``

    **Format :** int (natural number)

    **Description :** Number of temperature divisions.


- About replica 

    -  ``nprocs_per_replica``

    **Format :** int (natural number)

    **Description :** The number of processes for the replica. Default value = 1.

    -  ``nreplicas``

    **Format :** int (natural number)

    **Description :** The number of replicas.


- Others

    -  ``RXtrial_frequency``

    **Format :** int (natural number)

    **Description :** The number of Monte Carlo steps for replica exchange. Default = 1.


    -  ``sample_frequency``

    **Format :** int (natural number)

    **Description :**     The number of Monte Carlo steps for observation of physical quantities. Default value = 1.

    -  ``print_frequency``

    **Format :** int (natural number)

    **Description :**     The number of Monte Carlo steps for saving physical quantities. Default value = 1.

    -  ``reload``

    **Format :** bool ("true" or "false")

    **Description :**     Whether to restart the calculation from the last step finished last time. Default value = false.
