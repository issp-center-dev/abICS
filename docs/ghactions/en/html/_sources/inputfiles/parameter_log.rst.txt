.. highlight:: none

[log] section
-------------

This section specifies the log file name and the log level.


Input Format
^^^^^^^^^^^^
Keywords and their values are specified by a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Keywords
^^^^^^^^^^

   -  ``level``

      **Format :** str

      **Description :**
      Logging level. The following levels are available.
      
      - ``debug``
      - ``info``
      - ``warning``
      - ``error``

   - ``console``

      **Format :** str

      **Description :**
      console output mode.

         - ``default`` will examine if MPI environment is available or not.
         - ``mpi`` for parallel environment in which rank numbers are shown in error log.
         - ``serial`` for serial environment.
         - ``none`` suppresses console output.

   - ``console_level``

      **Format :** str

      **Description :**
      Logging level for console output.

   - ``logfile_path``

      **Format :** str

      **Description :**
      Path to the log file. If not specified, logs will be send only to console.
      The parent directories will be automatically created if they are not present.

   - ``logfile_mode``

      **Format :** str

      **Description :**
      MPI log type.

      - ``master`` will output logs to fiile only from rank=0.
      - ``collect`` will write messages from all ranks to one file.
      - ``workers`` will open one log file for each process designated by its rank.
      - ``serial`` will not consider parallel environment.

   - ``logfile_level``

      **Format :** str

      **Description :**
      Logging level for log file output.

   - ``logfile_rank``

      **Format :** int or list of int

      **Description :**
      MPI ranks from which logs are written to file.
      If not specified, all ranks are taken account of.
