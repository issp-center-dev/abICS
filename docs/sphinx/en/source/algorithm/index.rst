***************************
Sampling algorithm
***************************
abICS is designed for combining parallel extended ensemble methods with
arbitrary energy calculators. At present, only the replica exchange
Monte Carlo method is implemented.

Replica exchange Monte Carlo method
------------------------------------
A disadvantage of the widely-used Metropolis Monte Carlo algorithm is
that it tends to get stuck in local minima.
The replica exchange approach aims to overcome this problem by
considering multiple copies, or replicas, of the system under study.
The algorithm may be described roughly as follows
(see references below for more accurate descriptions).
Monte Carlo sampling is performed on each replica independently at
varying temperatures. At preset intervals, the temperatures are
exchanged according to a Metropolis criterion that essentially
assigns lower temperatures to replicas that happen to have lower
energies. This allows an efficient sampling of the global configuration
space using replicas at higher temperatures and accurate sampling of
the local energy landscape at lower temperatures.

- Overview of abICS

  - `S. Kasamatsu and O. Sugino, J. Phys. Condens. Matter, 31, 085901 (2019) <https://iopscience.iop.org/article/10.1088/1361-648X/aaf75c/meta>`_.

- About Exchange Monte Carlo method

  - `K. Hukushima and K. Nemoto, J. Phys. Soc. Japan, 65, 1604 (1996) <https://journals.jps.jp/doi/abs/10.1143/JPSJ.65.1604>`_.
  - `R. Swendsen and J. Wang, Phys. Rev. Lett. 57, 2607 (1986) <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.57.2607>`_.



