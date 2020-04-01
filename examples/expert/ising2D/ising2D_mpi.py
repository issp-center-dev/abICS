import copy
from abics.mc import *
from abics.mc_mpi import *
from ising2D import ising2D, ising2D_config, observer


if __name__ == "__main__":
    # Each MPI process has different random seeds after RX_MPI_init()
    comm, nreplicas, nprocs_per_replica = RX_MPI_init()
    J = -1.0
    model = ising2D(J)
    size = 10
    nspin = size * size
    eqsteps = 2 ** 16 * 10  # nspin*1000
    nsteps = 2 ** 16 * 40  # nspin*1000
    sample_frequency = 1
    RXtrial_frequency = 2 ** 2
    config = ising2D_config(size, size)
    config.prepare_random()
    configs = [copy.deepcopy(config) for i in range(nreplicas)]
    # When TemperatureRX_MPI.run() is called, each MPI rank works on
    # configs[rank]. configs[rank] is different amongst ranks
    # due to different random seeds (see above).

    kTs = np.linspace(0.01, 5.0, nreplicas)
    parallelCalc = TemperatureRX_MPI(comm, CanonicalMonteCarlo, model, configs, kTs)
    parallelCalc.run(eqsteps, RXtrial_frequency)
    # parallelCalc.reload()
    obs_mean = parallelCalc.run(
        nsteps, RXtrial_frequency, sample_frequency, observer=observer()
    )

    if comm.Get_rank() == 0:
        for i in range(len(kTs)):
            print(
                kTs[i],
                "\t",
                "\t".join(
                    [str(obs_mean[i][j] / nspin) for j in range(len(obs_mean[i]))]
                ),
            )
            # pickle.dump(obs_all, open("obs_all.pickle", "wb"))
