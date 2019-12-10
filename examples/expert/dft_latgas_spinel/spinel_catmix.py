import numpy as np
import copy
from mpi4py import MPI

from pymatgen import Lattice, Structure, Element, PeriodicSite
from pymatgen.io.vasp import Poscar, VaspInput
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
from py_mc.mc import CanonicalMonteCarlo, grid_1D, observer_base
from py_mc.mc_mpi import RX_MPI_init, TemperatureRX_MPI
from py_mc.applications.latgas_abinitio_interface.model_setup import (
    group,
    defect_sublattice,
    config,
    dft_latgas,
    g_r,
)
from py_mc.applications.latgas_abinitio_interface.run_vasp_mpi import (
    test_runner,
    vasp_runner,
)


kB = 8.6173e-5
comm, nreplicas, nprocs_per_replica = RX_MPI_init()

################## RXMC parameters ################################
# specify temperatures for each replica, number of steps, etc.
kTstart = 500.0
kTend = 1500.0
kTs = kB * np.linspace(kTstart, kTend, nreplicas)
# kTstep = 1.1
# kTs = kB*np.array([kTstart*kTstep**i for i in range(nreplicas)])

# eqsteps = 2000 # Number of steps for equilibration.

# Set Lreload to True when restarting
Lreload = False
# Lreload = True
nsteps = 1000  # Number of steps for sampling
RXtrial_frequency = 2
sample_frequency = 1
print_frequency = 1
# specify grid for calculating g(r)
# dr = 0.01
# maxr = 5
# grid = grid_1D(dr, dr, maxr)


################### model setup ###############################
# we first choose a "model" defining how to perform energy calculations and trial steps
# on the "configuration" defined below
baseinput = VaspInput.from_directory("baseinput")
energy_calculator = vasp_runner(
    base_input_dir="./baseinput",
    path_to_vasp="/home/i0009/i000900/src/vasp.5.3/vasp.spawnready.gamma",
    nprocs_per_vasp=nprocs_per_replica,
    comm=MPI.COMM_SELF,
    perturb=0.1,
)
# energy_calculator = test_runner()
model = dft_latgas(energy_calculator, save_history=False)
##############################################################


##############  defect sublattice setup #######################

# The POSCAR file contains primitive cell of MgAl2O4 spinel structure
# with 8 O, 4 Al, and 2 Mg
spinel_str = Structure.from_file("POSCAR")
spinel_str.make_supercell([2, 2, 2])

# The "base structure" contains ion sites with no disorder
# In our case, we consider perfectly ordered O sublattice
base_str = spinel_str.copy()
base_str.remove_species(["Al", "Mg"])

# We consider disorder in cation sublattice
cation_str = spinel_str.copy()
cation_str.remove_species(["O"])
cation_sites = cation_str.frac_coords

numMg = int(cation_str.composition["Mg"])
numAl = int(cation_str.composition["Al"])

# Define groups that we want to sample;
# we only consider single-atom groups here
# V = group("V", [], []) # for vacancy
Al = group("Al", ["Al"])
Mg = group("Mg", ["Mg"])
cations = [Al, Mg]

# Define defect sublattices using fractional coordinates of sites and atom groups that
# occupy those sites
defect_sublattices = [defect_sublattice(cation_sites, cations)]
num_defects = [{"Mg": numMg, "Al": numAl}]

# Finally, prepare configuration definition
spinel_config = config(base_str, defect_sublattices, num_defects, [1, 1, 1])
spinel_config.shuffle()
configs = []
for i in range(nreplicas):
    configs.append(copy.deepcopy(spinel_config))


################### Observer definition #####################
class observer_spinel(observer_base):
    def __init__(self, Asite_struct, Bspecie):
        super(observer_spinel, self).__init__()
        self.Asite_struct = Asite_struct
        self.Bspecie = Bspecie
        self.site_matcher = StructureMatcher(
            ltol=0.1,
            primitive_cell=False,
            allow_subset=True,
            comparator=FrameworkComparator(),
            ignored_species=["O"],
        )

    def DOI(self, calc_state):
        asites = self.site_matcher.get_mapping(
            calc_state.config.structure, self.Asite_struct
        )
        # print asites
        # print spinel_config.structure
        # print spinel_config.Asite_struct
        x = 0
        for i in asites:
            if calc_state.config.structure.species[i] == Element(self.Bspecie):
                x += 1
        x /= float(len(asites))
        return x

    def logfunc(self, calc_state):
        return calc_state.energy, self.DOI(calc_state)


Asite_struct = cation_str.copy()
Asite_struct.remove_species(["Al"])
myobserver = observer_spinel(Asite_struct, "Al")

################### RXMC calculation #########################
RXcalc = TemperatureRX_MPI(comm, CanonicalMonteCarlo, model, configs, kTs)
# obs = RXcalc.run(eqsteps, RXtrial_frequency, sample_frequency, observer=myobserver, subdirs=True)
if Lreload:
    RXcalc.reload()
obs = RXcalc.run(
    nsteps,
    RXtrial_frequency,
    sample_frequency,
    print_frequency,
    observer=myobserver,
    subdirs=True,
)


if comm.Get_rank() == 0:
    print(obs)
