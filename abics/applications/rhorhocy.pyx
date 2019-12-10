#from scipy.stats import norm
import numpy as np
cimport numpy as np
from scipy.fftpack import fftn, ifftn
#import copy
#from pymatgen import Lattice, Structure, Element, PeriodicSite
#from pymatgen.io.vasp.inputs import Poscar
from mpi4py import MPI


def rho_autocorr(np.ndarray[np.float64_t,ndim=2] x, np.ndarray[np.float64_t,ndim=1] rho, structure, double dr, double rmax, comm, int nsim):
    lattice = structure.lattice
    cdef int numgrid = len(x)
    assert len(rho) == numgrid
    cdef double distmax = np.amax(structure.lattice.abc)

    cdef np.ndarray[np.float64_t,ndim=1] X = np.arange(0., distmax, dr)
    
    cdef np.ndarray[np.float64_t,ndim=1] autocorr = np.zeros(len(X))

    cdef int size = comm.Get_size()
    cdef int rank = comm.Get_rank()
    if rank == 0:
        print(X)
    assert size < numgrid
    if  numgrid%size != 0 or numgrid%nsim !=0:
        raise ValueError('mesh size must be divisible by the communicator size and nsim')
        
    cdef int numgrid_loc = numgrid//size
    cdef int numgrid_sim = numgrid//nsim

    
    #dist_bin_histogram = np.zeros(len(X),dtype=int)
    #bins = range(len(X)+1)
    cdef int i,j
    cdef np.ndarray[np.float64_t,ndim=2] dist
    cdef np.ndarray[np.float64_t,ndim=2] rhoi_rhoj
    cdef np.ndarray[np.int_t,ndim=2] dist_bin
    for i in range(numgrid_sim):
        dist = lattice.get_all_distances(
            x[rank*numgrid_loc:(rank+1)*numgrid_loc],
            x[i*nsim:(i+1)*nsim]
        )
        dist_bin = np.around(dist/dr).astype(int)
        #dist_bin_histogram += np.histogram(dist_bin,bins)[0]
        rhoi_rhoj = np.outer(
            rho[rank*numgrid_loc:(rank+1)*numgrid_loc],
            rho[i*nsim:(i+1)*nsim]
        )
        for j in range(int(rmax//dr)):
            autocorr[j] += np.sum((dist_bin==j)*rhoi_rhoj)

    comm.Allreduce(MPI.IN_PLACE, autocorr, op=MPI.SUM)
    #comm.Allreduce(MPI.IN_PLACE, dist_bin_histogram, op=MPI.SUM)
    
    X[0] = 1 # to avoid division by zero
    #autocorr*= lattice.volume/ \
    #             (4.0*np.pi*X*X*dr*numgrid*(numgrid-1))
    #dist_bin_histogram[0] = 1
    #autocorr /= dist_bin_histogram
    autocorr[0] = 0
    return autocorr[0:int(rmax//dr)]
    
    

def rho_autocorr_fft(mesh,x, rho, structure, dr, rmax): #, comm, nsim):
    numgrid = len(x)
    lattice = structure.lattice
    dist = lattice.get_all_distances([0.,0.,0.],x)
    dist = np.reshape(dist, (mesh[0],mesh[1],mesh[2]),order='F')
    dist_bin = np.around(dist/dr).astype(int)

    cdef int i,j,k, mesh0, mesh1, mesh2
    cdef np.ndarray[np.complex128_t,ndim=3] c_q, rho_K
    rhoxyz = np.reshape(rho,(mesh[0],mesh[1],mesh[2]),order='F')
    rho_K = fftn(rhoxyz)
    mesh0 = mesh[0]
    mesh1 = mesh[1]
    mesh2 = mesh[2]
    c_q = np.zeros((mesh[0],mesh[1],mesh[2]),dtype=complex)
    for i in range(mesh0):
        for j in range(mesh1):
            for k in range(mesh2):
                c_q[i,j,k] = rho_K[i,j,k] * rho_K[-i,-j,-k]
    c_r = ifftn(c_q)
    
    distmax = np.amax(structure.lattice.abc)
    X = np.arange(0., distmax, dr)
    autocorr = np.zeros(len(X))
    X[0] = 1
    for j in range(int(rmax//dr)):
        autocorr[j] += np.sum((dist_bin==j)*c_r.real)
    autocorr*= lattice.volume/ \
                 (4.0*np.pi*X*X*dr*numgrid*(numgrid-1))
    autocorr[0] = 0
 
    
    return autocorr[0:int(rmax//dr)]
