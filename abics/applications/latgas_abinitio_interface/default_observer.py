from abics.mc import observer_base


class default_observer(observer_base):
    """
    Default observer.

    Attributes
    ----------
    minE : float
        Minimum of energy
    """
    def __init__(self, comm, Lreload=False):
        """

        Parameters
        ----------
        comm: mpi4py.MPI.Intracomm
            MPI communicator
        Lreload: bool
            Reload or not
        """
        super(default_observer, self).__init__()
        self.minE = 100000.0
        myrank = comm.Get_rank()
        if Lreload:
            minEfi = open(str(myrank) + "/minEfi.dat", "r")
            self.minE = float(minEfi.readlines()[-1])
            minEfi.close()
            obs_fi = open(str(myrank) + "/obs.dat", "r")
            self.lprintcount = int(obs_fi.readlines()[-1].split()[0]) + 1

    def logfunc(self, calc_state):
        """

        Parameters
        ----------
        calc_state: MCalgo
        Object of Monte Carlo algorithm
        Returns
        -------
        calc_state.energy : float
        Minimum energy
        """
        if calc_state.energy < self.minE:
            self.minE = calc_state.energy
            minEfi = open("minEfi.dat", "a")
            minEfi.write(str(self.minE) + "\n")
            calc_state.config.structure.to(fmt="POSCAR", filename="minE.vasp")
        return calc_state.energy

    def writefile(self, calc_state):
        """

        Parameters
        ----------
        calc_state: MCalgo
        Object of Monte Carlo algorithm

        """
        calc_state.config.structure.to(
            fmt="POSCAR", filename="structure." + str(self.lprintcount) + ".vasp"
        )
