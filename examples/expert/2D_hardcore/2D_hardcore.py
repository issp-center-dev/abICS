import numpy as np
import random as rand
import sys

from py_mc.mc import CanonicalMonteCarlo, grid_1D
from model_setup import *


if __name__ == "__main__":
    L = 1.0
    rdisc = 1.0 / 14.0
    Ndisc = 20
    dr = 0.01
    maxstep = rdisc / 2
    kT = 1

    grid = grid_1D(dr, dr, L / 2.0)
    config = HC_2D_config(L, Ndisc, rdisc)
    model = HC_2D(maxstep)
    config.prepare_ordered()
    print(config.coors)
    print(g_r(config, grid))
    calc = CanonicalMonteCarlo(model, kT, config, grid=grid)
    plot_fig(calc, 10000)
    # calc.run(10000, observefunc=observables)
    # print(config.coors)

    """calc = Canonical
    
    J = -1.0
    kT = abs(J) * 1.0
    size = 10
    eqsteps = 100000
    mcsteps = 1000000
    sample_frequency = size*size
    config = ising2D_config(size,size)
    config.prepare_random()
    model = ising2D(J)

    for kT in [5]: #np.arange(5, 0.5, -0.05):
        energy_expect = 0
        magnet_expect = 0
        
        kT = abs(J)*kT        

        #print config        
        calc = CanonicalMonteCarlo(model, kT, config)
        calc.run(eqsteps)

        mcloop = mcsteps//sample_frequency
        for i in range(mcloop):
            calc.run(sample_frequency)
            #print model.energy(config), model.magnetization(config)
            current_config = calc.config
            energy_expect += model.energy(current_config)
            magnet_expect += abs(model.magnetization(current_config))
        print(kT, energy_expect/mcloop, magnet_expect/mcloop)
        sys.stdout.flush()
    #calc.run(100000)
    #print config
    
        
"""
