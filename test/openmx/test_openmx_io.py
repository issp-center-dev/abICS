import subprocess
import sys
import os

if __name__ == '__main__':
    sys.path.append('../../')
    from abics.applications.latgas_abinitio_interface.openmx import OpenMXSolver    
    #Set solver
    path_to_solver = ""
    solver = OpenMXSolver(path_to_solver)
    print ("solver name is {}. ".format(solver.name()))

    base_dir = "/Users/k-yoshimi/Dropbox/PycharmProjects/abICS/test/openmx/sample"
    #Read base file
    print("Read base file.")
    input = solver.input
    input.from_directory(base_dir)
    #Execute OpenMX
    print("Write input file.")
    input_dir = "test"
    input.write_input(input_dir)
    print("Execute OpenMX.")
    cmd = "openmx {}.dat".format(os.path.join(input_dir, input.base_openmx_input["System.Name"][0]))
    subprocess.call(cmd.split())

    #Read Output file
    print("Read Output file.")
    output = solver.output
    phys = output.get_results("./")
    print("Update input file.")
    input.update_info_by_structure(phys.structure)
    print("Write input file.")
    input_dir = "test1"
    input.write_input(input_dir)
    print("Execute OpenMX.")
    cmd = "openmx {}.dat".format(os.path.join(input_dir, input.base_openmx_input["System.Name"][0]))
    subprocess.call(cmd.split())
    # Check the results of test/met.dat# and test1/met.dat# are almost same.

    #Check seldyn_arr
    phys.structure.add_site_property("seldyn", [[False, False, False]] * len(phys.structure.sites))
    input.update_info_by_structure(phys.structure)
    print("Write input file.")
    input_dir = "test2"
    input.write_input(input_dir)
    print("Execute OpenMX.")
    cmd = "openmx {}.dat".format(os.path.join(input_dir, input.base_openmx_input["System.Name"][0]))
    subprocess.call(cmd.split())
    # Check the coordinates of test1/met.dat# and test2/met.dat# are same.

    # test subprocess mode
    from abics.applications.latgas_abinitio_interface.run_base_mpi import runner
    from mpi4py import MPI
    nprocs_per_replica = 1
    energy_calculator = runner(
        base_input_dir=base_dir,
        Solver=solver,
        nprocs_per_solver=nprocs_per_replica,
        comm=MPI.COMM_SELF,
        solver_run_scheme = "subprocess"
    )
    energy_calculator.submit(structure=phys.structure, output_dir="./")
