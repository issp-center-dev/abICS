.. highlight:: none

[observer] section
-------------------------------

Specify the type of physical quantity to be acquired.
The example is shown as follows:

  ::
     
    [observer]
    type = 'default'

Input Format
^^^^^^^^^^^^^
Specify a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Key words
^^^^^^^^^^

 -  ``type``

    **Format :** str 

    **Description :**
    Specify a physical quantity set.
    
    - "default"
        
        - Default setting. Get the energy and the atom group in each coordinate.
