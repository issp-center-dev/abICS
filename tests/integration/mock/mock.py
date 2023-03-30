def calc_energy(st):
    ene = 0.0
    st_local = st.copy()
    dm = st_local.distance_matrix
    n = len(st_local)
    an_mean = 0
    for i in range(n):
        an_mean += st_local.species[i].number
    an_mean /= n
    for i in range(n):
        an_i = st_local.species[i].number - an_mean
        for j in range(i + 1, n):
            an_j = st_local.species[j].number - an_mean
            ene += (an_i * an_j) / (dm[i, j] ** 2)
    return ene
