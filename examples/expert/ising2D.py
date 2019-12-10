import numpy as np
import random as rand
import sys

from py_mc.mc import model, CanonicalMonteCarlo, binning, observer_base


class ising2D(model):
    """This class defines the 2D ising model"""

    model_name = "ising2D"

    def __init__(self, J):
        self.J = J

    def energy(self, ising2D_config):
        """ Calculate total energy of the 2D Ising model"""
        e = 0.0
        config = ising2D_config.config  # config should be 2D numpy array
        for i in range(ising2D_config.lenX):
            for j in range(ising2D_config.lenY):
                e += self.J * config[i, j] * (config[i - 1, j] + config[i, j - 1])
        return e

    def magnetization(self, ising2D_config):
        """Calculate magnetization"""
        return ising2D_config.config.sum()

    def trialstep(self, ising2D_config, energy):
        # energy is just a placeholder and isn't used
        # here
        # choose x,y randomly
        x = rand.randrange(ising2D_config.lenX)
        y = rand.randrange(ising2D_config.lenY)
        # Position of flipped spin to be used by newconfig():
        dconfig = [x, y]
        config = ising2D_config.config

        # Calculate energy change if spin flips at x,y
        x1 = x + 1
        if x == ising2D_config.lenX - 1:
            x1 = 0
        y1 = y + 1
        if y == ising2D_config.lenY - 1:
            y1 = 0
        dE = (
            -2.0
            * self.J
            * config[x, y]
            * (config[x - 1, y] + config[x1, y] + config[x, y - 1] + config[x, y1])
        )
        # print dconfig, dE
        return dconfig, dE

    def newconfig(self, ising2D_config, dconfig):
        """Construct the new configuration after the trial step is accepted"""
        ising2D_config.config[dconfig[0], dconfig[1]] *= -1
        return ising2D_config


class ising2D_config:
    """This class defines the 2D Ising model configuration"""

    def __init__(self, lenX, lenY):
        self.lenX = lenX
        self.lenY = lenY
        self.config = np.empty([lenX, lenY])

    def prepare_random(self):
        for i in range(self.lenX):
            for j in range(self.lenY):
                if rand.random() >= 0.5:
                    self.config[i, j] = 1
                else:
                    self.config[i, j] = -1

    def __str__(self):
        s = ""
        for i in range(self.lenX):
            for j in range(self.lenY):
                if self.config[i, j] < 0:
                    s += "-"
                else:
                    s += "+"
            s += "\n"
        return s


class observer(observer_base):
    def __init__(self):
        super(observer, self).__init__()
        self.energy_obs = []
        self.magnet_obs = []

    def logfunc(self, calc_state):
        energy = calc_state.energy
        magnetization = calc_state.model.magnetization(calc_state.config)
        absmag = abs(magnetization)
        # Save observations at each step
        self.energy_obs.append(energy)
        self.magnet_obs.append(absmag)
        return energy, absmag


if __name__ == "__main__":
    J = -1.0
    kT = abs(J) * 5.0
    size = 10
    nspin = size * size
    nstep_param = 12
    eqsteps = 2 ** nstep_param * 10  # equilibration steps
    mcsteps = 2 ** nstep_param * 40  # measurement steps
    sample_frequency = 1  # we sample every step
    print_frequency = (
        100000000
    )  # the frequency for printing observer.logfunc() to obs.dat

    binlevel = int(np.log2(mcsteps // 30))
    config = ising2D_config(size, size)
    config.prepare_random()
    model = ising2D(J)
    binning_fileE = open("binningE.dat", "w")
    binning_fileM = open("binningM.dat", "w")
    for kT in np.linspace(5.0, 0.01, 10):
        kT = abs(J) * kT
        calc = CanonicalMonteCarlo(model, kT, config)
        calc.run(eqsteps)  # run equilibration steps

        myobserver = observer()  # provide a new observer instance
        obs = calc.run(
            mcsteps, sample_frequency, print_frequency, myobserver  # Run sampling steps
        )

        # binning analysis for putting error bars on data
        error_estimateE = binning(np.asarray(myobserver.energy_obs) / nspin, binlevel)
        error_estimateM = binning(np.asarray(myobserver.magnet_obs) / nspin, binlevel)
        binning_fileE.write("\n".join([str(x) for x in error_estimateE]) + "\n\n\n")
        binning_fileM.write("\n".join([str(x) for x in error_estimateM]) + "\n\n\n")
        print(
            kT,
            "\t",
            "\t".join([str(x / nspin) for x in obs]),
            np.max(error_estimateE),
            np.max(error_estimateM),
        )
        sys.stdout.flush()
        # model = calc.model
        # config = calc.config
        # print(config)
