import sys
import os

from pymatgen.core import Structure

def calc_energy(st: Structure) -> float:
    ene = 0.0
    st_local = st.copy()
    dm = st_local.distance_matrix
    n = len(st_local)
    an_mean = 0.0
    for i in range(n):
        an_mean += st_local.species[i].number
    an_mean /= n
    for i in range(n):
        an_i = st_local.species[i].number - an_mean
        for j in range(i + 1, n):
            an_j = st_local.species[j].number - an_mean
            ene += (an_i * an_j) / (dm[i, j] ** 2)
    return ene

if __name__ == "__main__":
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    stfile = os.path.join(output_dir, "structure.vasp")
    st = Structure.from_file(stfile)
    energy = calc_energy(st)
    with open(os.path.join(output_dir, "energy.dat"), "w") as f:
        f.write("{:.15f}\n".format(energy))
