import os, sys
import numpy as np
import numpy.random as random
import scipy.constants as constants
import copy
import toml

from abics.applications.latgas_abinitio_interface.defect import (
    DFTConfigParams,
    defect_config,
)
from abics.applications.latgas_abinitio_interface.model_setup import DFTLatticeGas
from abics.sampling.mc import WeightedCanonicalMonteCarlo
from abics.observer import ObserverBase

import logging

logger = logging.getLogger("main")


class EnergyCalculator:
    def __init__(self, params):
        self.J = params["J"]
        self.mu = params["mu"]
        self.Lx = params["Lx"]
        self.Ly = params["Ly"]
        self.kbT = params["kbT"]

    def calc(self, config):
        # make backup (no-relaxed state)
        config.structure_norel = config.structure.copy()

        # interprete latgas_rep
        state = np.array(
            [(1 if v[0] == "A" else 0) for v in config.defect_sublattices[0].latgas_rep]
        ).reshape(self.Lx, self.Ly)
        # E = J Sum_<ij> e_i e_j - mu N  for e_i = {0, 1}
        state1 = (np.roll(state, 1, axis=0) + np.roll(state, 1, axis=1)) * state

        # e0 = 0.0 + np.sum(state1) * self.J - np.sum(state) * self.mu
        e0 = 0.0 + np.sum(state1) * self.J

        logger.debug("EnergyCalculator.calc: e = {}".format(e0))
        return e0


class Observer(ObserverBase):
    def __init__(self):
        super().__init__()
        self.energy_obs = []
        self.density_obs = []

    def logfunc(self, calc_state):
        energy = calc_state.energy
        num = calc_state.model.density(calc_state.config)

        self.energy_obs.append(energy)
        self.density_obs.append(num)

        logger.debug("Observer.logfunc: energy={}, density={}".format(energy, num))
        return energy, num


class LatticeGas(DFTLatticeGas):
    def __init__(
        self,
        abinitio_run,
        enable_grandcanonical=True,
        gc_ratio=1.0,
    ):
        super().__init__(
            abinitio_run=abinitio_run,
            save_history=False,
            l_update_basestruct=False,
            check_ion_move=False,
            ion_move_tol=0.7,
            enable_grandcanonical=enable_grandcanonical,
            gc_ratio=gc_ratio,
            debug=False,
        )
        logger.debug(">>> LatticeGas()")

    def energy(self, config):
        logger.debug(">>> LatticeGas.energy()")
        if config.energy is None:
            config.energy = self.internal_energy(config)
            if self.enable_grandcanonical:
                config.energy -= self._calc_muN_term(config)
        return config.energy

    def internal_energy(self, config):
        logger.debug(">>> LatticeGas.internal_energy()")
        if config.internal_energy is None:
            config.internal_energy = self._calc_energy(config)
        return config.internal_energy

    def _calc_energy(self, config):
        logger.debug(">>> LatticeGas._internal_energy()")
        ev = self.abinitio_run.calc(config)
        return ev

    def density(self, config):
        state = np.array(
            [(1 if v[0] == "A" else 0) for v in config.defect_sublattices[0].latgas_rep]
        )
        num = np.sum(state)
        return num


def main(params):
    params_config = params.get("config", {})
    params_sampling = params.get("sampling", {})

    # physical constants
    kb = constants.Boltzmann
    m_ = 1.0e-3 / constants.Avogadro
    h_ = constants.Planck
    eV = constants.elementary_charge

    # parameters
    T = params_config.get("T", 300)
    kbT = kb * T
    mu0 = -kbT * np.log((2.0 * np.pi * m_ * kbT / h_**2) ** (3 / 2) * kbT)
    Eads = params_config.get("Eads", -0.4) * eV
    J = params_config.get("J", 0.0) * eV

    Lx, Ly = params_config.get("L", [5, 5])
    nspin = Lx * Ly
    eqsteps = params_sampling.get("thermalization_steps", 1) * nspin
    mcsteps = params_sampling.get("sample_steps", 1) * nspin
    sample_frequency = params_sampling.get("sample_frequency", 1)
    print_frequency = params_sampling.get("print_frequency", 1)

    pstart = params_sampling.get("p_start", 1.0e4)
    pend = params_sampling.get("p_end", 1.0e-4)
    psteps = params_sampling.get("p_steps", 1)
    pscale = params_sampling.get("p_scale", "linear")

    if pscale == "log" or pscale == "logarithm":
        ptable = np.logspace(np.log10(pstart), np.log10(pend), psteps, endpoint=True)
    else:
        ptable = np.linspace(pstart, pend, psteps, endpoint=True)

    # initialize random number sequence
    random_seed = params_sampling.get("random_seed", 123456789)
    random.seed(random_seed)

    conf_start = params_sampling.get("conf_start", 0.0)
    if conf_start < 0:
        nstart = random.randint(nspin + 1)
    else:
        nstart = round(nspin * conf_start)
        if nstart < 0:
            nstart = 0
        elif nstart > nspin:
            nstart = nspin

    conf_base = {
        "unitcell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "base_structure": [
            {
                "type": "O",
                "coords": [[0, 0, 0]],
            },
        ],
        "defect_structure": [
            {
                "coords": [
                    [1.0 * i / (Lx - 1), 1.0 * j / (Ly - 1), 0.0]
                    for i in range(Lx)
                    for j in range(Ly)
                ],
                "groups": [
                    {
                        "name": "A",
                        "num": 0,
                    }
                ],
            },
        ],
        "chemical_potential": [
            {
                "species": ["A"],
                "mu": 0.0,
            },
        ],
        "grandcanonical_move": [
            {"species": "A"},
        ],
    }

    # # initialize random number sequence
    # random.seed(random_seed)

    # --------
    logger.info("parameters:")
    # ______________123456789012345678901234567890
    logger.info("  T                        = {}".format(T))
    logger.info("  Eads                     = {}".format(Eads / eV))
    logger.info("  J                        = {}".format(J / eV))
    logger.info("  Lx,Ly                    = {}".format([Lx, Ly]))
    logger.info("  thermalization_steps     = {}".format(eqsteps // nspin))
    logger.info("  sample_steps             = {}".format(mcsteps // nspin))
    logger.info("  sample_frequency         = {}".format(sample_frequency))
    logger.info("  print_frequency          = {}".format(print_frequency))
    logger.info("  random_seed              = {}".format(random_seed))
    logger.info("  initial occupancy        = {}".format(nstart / nspin))

    # logger.info("ptable = {}".format(ptable))

    # --------
    for p in ptable:
        mu = mu0 + kbT * np.log(p)

        logger.info("step p={}, mu={}".format(p, mu))
        logger.info("step p={}, mu-Eads={}".format(p, mu - Eads))

        # create parameter set
        params_config_input = copy.deepcopy(conf_base)
        params_config_input["defect_structure"][0]["groups"][0]["num"] = nstart
        params_config_input["chemical_potential"][0]["mu"] = mu - Eads

        params_config = DFTConfigParams.from_dict(params_config_input)

        # create config
        config = defect_config(params_config)

        calc_params = {"Lx": Lx, "Ly": Ly, "J": J, "mu": mu - Eads, "kbT": kbT}
        energy_calculator = EnergyCalculator(calc_params)

        model = LatticeGas(energy_calculator, enable_grandcanonical=True, gc_ratio=1.0)

        mc = WeightedCanonicalMonteCarlo(model, kbT, config)

        # thermalization loop
        mc.run(eqsteps)

        # observer
        myobserver = Observer()

        # main loop
        obs = mc.run(mcsteps, sample_frequency, observer=myobserver)

        print(p, *([str(x / nspin) for x in obs]), sep="\t")

        # error_estimate = binning(myobserver.density_obs, 10)

        # model = mc.model
        # config = mc.config


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--verbose", help="increase verbosity", action="count", default=0
    )
    parser.add_argument(
        "-q", "--quiet", help="decrease verbosity", action="count", default=0
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default="input.toml",
        help="parameter file in TOML format",
    )
    args = parser.parse_args()

    log_level = logging.INFO - (args.verbose - args.quiet) * 10
    logging.basicConfig(level=log_level)

    try:
        params = toml.load(args.input_file)
    except Exception as e:
        print(e)
        sys.exit(1)

    main(params)
