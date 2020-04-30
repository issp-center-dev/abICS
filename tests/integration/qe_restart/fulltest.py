import os
import shutil
import sys
import unittest

from mpi4py import MPI
import numpy as np

sys.path.insert(0, "../../../")

from abics.scripts.main import main_impl


class TestAbics(unittest.TestCase):
    def test_hoge(self):
        if MPI.COMM_WORLD.Get_rank() == 0:
            for i in (0, 1):
                if os.path.exists(str(i)):
                    shutil.rmtree(str(i))
                if os.path.exists("{}_4step".format(i)):
                    shutil.rmtree("{}_4step".format(i))
        MPI.COMM_WORLD.Barrier()

        main_impl("input_4steps.toml")
        if MPI.COMM_WORLD.Get_rank() == 0:
            ref_T = np.zeros((2, 4))
            ref_E = np.zeros((2, 4))
            for i in (0, 1):
                with open(os.path.join(str(i), "obs.dat")) as f:
                    for j, line in enumerate(f):
                        words = line.split()
                        ref_T[i, j] = float(words[1])
                        ref_E[i, j] = float(words[2])
                        if j == 3:
                            break
                shutil.move(str(i), "{}_4step".format(i))

        main_impl("input_2steps.toml")
        main_impl("input_2steps_restart.toml")

        if MPI.COMM_WORLD.Get_rank() == 0:
            res_T = np.zeros((2, 4))
            res_E = np.zeros((2, 4))
            for i in (0, 1):
                with open(os.path.join(str(i), "obs.dat")) as f:
                    for j, line in enumerate(f):
                        words = line.split()
                        res_T[i, j] = float(words[1])
                        res_E[i, j] = float(words[2])
                        if j == 3:
                            break
            self.assertTrue(np.allclose(res_T, ref_T))
            self.assertTrue(np.allclose(res_E, ref_E))


if __name__ == "__main__":
    unittest.main()
