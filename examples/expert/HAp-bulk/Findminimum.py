import numpy as np
import random as rand
import sys, os
import copy
import pickle


if __name__ == "__main__":
    energy_lst = pickle.load(open("energy_reps.pickle", "rb"))
    min_id = np.argmin(energy_lst)
    print(min_id)
    reps = pickle.load(open("latgas_reps.pickle", "rb"))
    min_rep = reps[min_id]
    print(min_rep)
