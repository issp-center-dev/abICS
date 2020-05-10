import os
import unittest
import sys

from mpi4py import MPI
import numpy as np

sys.path.insert(0, "../../../")

from abics.scripts.main import main_impl


class TestAbics(unittest.TestCase):
    def test_hoge(self):
        main_impl("input.toml")
        if MPI.COMM_WORLD.Get_rank() == 0:
            ref_T = np.array([[0.103407999144, 0.103407999144], [0.08617333262, 0.08617333262]])
            ref_E = np.array(
                [[-5963.873264468331, -5963.8732643856365], [-5963.871843725528, -5963.873264468331]]
            )
            T = np.zeros((2, 2))
            E = np.zeros((2, 2))
            for i in (0, 1):
                with open(os.path.join(str(i), "obs.dat")) as f:
                    for j, line in enumerate(f):
                        words = line.split()
                        T[i, j] = float(words[1])
                        E[i, j] = float(words[2])
                        if j == 1:
                            break

            # print()
            # print(T)
            # print(E)
            # print()
            # print(ref_T)
            # print(ref_E)
            # print()

            self.assertTrue(np.allclose(T, ref_T))
            self.assertTrue(np.allclose(E, ref_E))


if __name__ == "__main__":
    unittest.main()
