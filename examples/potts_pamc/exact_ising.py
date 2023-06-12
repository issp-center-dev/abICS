# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

"""
Exact solution of the internal energy U,
the Helmholtz free energy logZ, and the squared magnetization m2 of
Q=2 Potts model on square lattice (equivalent to the Ising model).
Equations are taken from Wikipedia (2022-11-07)
https://en.wikipedia.org/wiki/Square_lattice_Ising_model
"""

import numpy as np
from scipy.integrate import quad


def Ising_m2(beta: float, J: float = 1.0) -> float:
    if beta == 0.0:
        return 0.0
    res = 1.0 - np.sinh(2.0 * beta * J) ** (-4)
    if res < 0.0:
        return 0.0
    else:
        return res**0.25


def Ising_U(beta: float, J: float = 1.0) -> float:
    assert beta >= 0.0
    if beta == 0.0:
        return 0.0
    K = beta * J
    k: float = 1.0 / (np.sinh(2 * K) ** 2)

    def integrant(t: float) -> float:
        return 1.0 / np.sqrt(1.0 - 4 * k / ((1 + k) * (1 + k)) * np.sin(t) ** 2)

    res: float = quad(integrant, a=0.0, b=0.5 * np.pi)[0]
    res *= (2.0 / np.pi) * ((2 * np.tanh(2 * K) ** 2) - 1)
    res += 1.0
    res /= -np.tanh(2 * K)
    return J * res


def Ising_mbF(beta: float, J: float = 1.0) -> float:
    assert beta >= 0.0
    if beta == 0.0:
        return np.log(2.0)
    K = beta * J
    L = K
    k: float = 1.0 / (np.sinh(2 * K) * np.sinh(2 * L))

    def integrant(t: float) -> float:
        res = np.cosh(2 * K) * np.cosh(2 * L)
        res += np.sqrt(1.0 + k * k - 2 * k * np.cos(2 * t)) / k
        return np.log(res)

    res = quad(integrant, a=0.0, b=np.pi)[0] / (2 * np.pi)
    return res + 0.5 * np.log(2.0)


def Potts_m2(beta: float, J: float = 1.0) -> float:
    return 0.25 * Ising_m2(beta, J=0.5 * J)


def Potts_U(beta: float, J: float = 1.0) -> float:
    return Ising_U(beta, J=0.5 * J) - J


def Potts_mbF(beta: float, J: float = 1.0) -> float:
    return Ising_mbF(beta, J=0.5 * J) + beta * J


if __name__ == "__main__":
    # Ts = np.concatenate((np.linspace(0.5, 2.0, num=101), np.array([1e8])))
    betas = np.linspace(0.0, 2.0, num=201)
    print("# beta U logZ m2")
    for beta in betas:
        u = Potts_U(beta, J=1.0)
        mbf = Potts_mbF(beta, J=1.0)
        m2 = Potts_m2(beta, J=1.0)
        print(beta, u, mbf, m2)
