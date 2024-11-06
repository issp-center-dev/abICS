.. highlight:: none

[sampling] section
-------------------------------

Specify the parameters of the Monte Carlo (MC) sampling method, such as the number of replicas, the temperature range, and the number of Monte Carlo steps.
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

- About sampling method

   - ``sampler``

       **Format :** string

       **Description :**
       Relica exchange MC method ("RXMC") or population annealing MC method ("PAMC").

- About temperatures

   - Specify temperature points by using ``kTs`` or ``kTstart``, ``kTend``, and ``kTnum`` (lineary spaced).
     If ``kTs`` is specified, the others will be ignored.
   - Temperatures should be given in the unit of Kelvin.

   - ``kTs``

       **Format :** list of float (>0)

       **Description :**
       Temperature points.
       When ``sampler = "RXMC"``, the number of temperature points should equal to ``nreplicas``.

   - ``kTstart``

       **Format :** float (>0)

       **Description :**
       Minimum temperature.

   - ``kTend``

       **Format :** float (>0)

       **Description :**
       Maximum temperature.

   - ``kTnum`` (Only for PAMC)

       **Format :** int (>0)

       **Description :**
       The number of temperature points.
       When ``sampler = "RXMC"``, the number of temperature points will equal to ``nreplicas``.

    - ``linspace_in_beta``

       **Format :** true or false

       **Description :**
       If true, temperature points are generated in the inverse temperature space with equal intervals.
       If false, temperature points are generated in the temperature space with equal intervals.
       Default value = false.

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

   - ``nsteps_between_annealing`` (Only for ``sampler = "PAMC"``)

       **Format :** int (natural number)

       **Description :** Number of Monte Carlo steps for each temperature.
  
   - ``RXtrial_frequency`` (Only for ``sampler = "RXMC"``)

       **Format :** int (natural number)

       **Description :** The interval for performing replica exchange trials. For example, setting this value to 1 means that replica exchange is attempted at every Monte Carlo step, while setting this to 2 means that exchange is attempted at every second step. Default = 1.

   - ``resample_frequency`` (Only for ``sampler = "PAMC"``)

       **Format :** int (natural number)

       **Description :** The interval for performing replica resampling. For example, setting this value to 1 means that replica resampling is attempted after every temperature lowering, while setting this to 2 means that resampling is attempted at every second step. Default = 1.

   - ``sample_frequency``

       **Format :** int (natural number)

       **Description :** The interval for observation of physical quantities. Default value = 1.

   - ``print_frequency``

       **Format :** int (natural number)

       **Description :** The interval for saving physical quantities. Default value = 1.

   - ``reload``

       **Format :** bool ("true" or "false")

       **Description :** Whether to restart a prior calculation from the last step finished last time. Default value = false.

    - ``throw_out``

       **Format :** int or float

       **Description :** The number (int) or ratio (float) of measurements to be thrown out as thermalization in the process of the evaluation of expectation values. Default value = 0.5 .

    - ``enable_grandcanonical``

       **Format :** bool ("true" or "false")

       **Description :** Whether to allow grand canonical sampling. Default value = false.

    - ``gc_ratio``

       **Format :** float

       **Description :** The ratio of the grand canonical update that changes the number of elements among the trials of configuration updates, when the grand canonical sampling is turned on. Default value = 0.3 .
