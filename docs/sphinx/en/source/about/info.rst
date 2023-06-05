.. pyMC documentation master file, created by
   sphinx-quickstart on Wed Jul 31 13:13:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

What is abICS ?
------------------------------------------
abICS is a software framework for training a machine learning model to
reproduce first-principles energies and then using the model to perform
configurational sampling in disordered systems.
Specific emphasis is placed on multi-component solid state systems such as metal and oxide alloys.
The current version of abics can use neural network models implemented in aenet to be used as 
the machine learning model. As of this moment, abICS can also generate Quantum Espresso, VASP, 
and OpenMX input files for obtaining the reference training data for the machine learning model.
For the sampling algorithms, abICS implements the extended Monte Carlo methods, namely, the replica exchange Monte Carlo method (RXMC) and the population annealing Monte Carlo method (PAMC).
In addition, as a beta version, the grand canonical sampling is supported.


Developers
------------------------------------------
abICS is developed by the following members.

- ver. 2.0
   - Shusuke Kasamatsu (Yamagata University)
   - Yuichi Motoyama (Institute for Solid State Physics, Univ. of Tokyo)
   - Kazuyoshi Yoshimi (Institute for Solid State Physics, Univ. of Tokyo)
   - Tatsumi Aoyama (Institute for Solid State Physics, Univ. of Tokyo)
   - Osamu Sugino (Institute for Solid State Physics, Univ. of Tokyo)

- ver. 1.0
   - Shusuke Kasamatsu (Yamagata University)
   - Yuichi Motoyama (Institute for Solid State Physics, Univ. of Tokyo)
   - Kazuyoshi Yoshimi (Institute for Solid State Physics, Univ. of Tokyo)
   - Yoshiyuki Yamamoto (Institute for Solid State Physics, Univ. of Tokyo)
   - Osamu Sugino (Institute for Solid State Physics, Univ. of Tokyo)
   - Taisuke Ozaki (Institute for Solid State Physics, Univ. of Tokyo)
   
Version information
------------------------------------------

- ver. 2.1.0    : 2023/06/08.
- ver. 2.0.1    : 2022/11/04.
- ver. 2.0      : 2022/06/24.
- ver. 1.0      : 2020/05/01.
- ver. 1.0-beta : 2020/03/31.
- ver. 0.1      : 2019/12/10.


License
--------------

This package is distributed under GNU General Public License version 3 (GPL v3) or later.


Copyright
--------------

*(c) 2019- The University of Tokyo. All rights reserved.*

This software was developed with the support of \"*Project for advancement of software usability in materials science*\" of The Institute for Solid State Physics, The University of Tokyo. 
     
