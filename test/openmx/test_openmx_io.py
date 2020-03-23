import subprocess
import sys
import os



if __name__ == '__main__':
    sys.path.append('../../')
    from abics.applications.latgas_abinitio_interface.openmx import OpenMXSolver    
    #Set solver
    # path_to_solver = "openmx"
    path_to_solver = os.path.expanduser("~/build/openmx3.9/work/openmx")
    solver = OpenMXSolver(path_to_solver)
    print ("solver name is {}. ".format(solver.name()))

    #base_dir = "/Users/k-yoshimi/Dropbox/PycharmProjects/abICS/test/openmx/sample"
    base_dir =  os.path.join(os.path.abspath(os.path.dirname(__file__)), "sample")
    #Read base file
    print("Read base file.")
    input = solver.input
    input.from_directory(base_dir)
    #Execute OpenMX
    print("Write input file.")
    input_dir = "test"
    input.write_input(input_dir)
    print("Execute OpenMX.")
    cmd = "{} {}.dat".format(path_to_solver, os.path.join(input_dir, input.base_openmx_input["System.Name"][0]))
    # subprocess.call(cmd.split())

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
    cmd = "{} {}.dat".format(path_to_solver, os.path.join(input_dir, input.base_openmx_input["System.Name"][0]))
    # subprocess.call(cmd.split())
    #Check the results of test/met.dat# and test1/met.dat# are almost same.

    #Check seldyn_arr
    print("Add seldyn_arr and check the result using subprocess scheme")
    phys.structure.add_site_property("seldyn", [[False, False, False]] * len(phys.structure.sites))
    input.update_info_by_structure(phys.structure)
    # test subprocess mode and seldyn
    from abics.applications.latgas_abinitio_interface.run_base_mpi import runner
    from mpi4py import MPI
    nprocs_per_replica = 1
    energy_calculator = runner(
        base_input_dir=base_dir,
        Solver=solver,
        nprocs_per_solver=nprocs_per_replica,
        comm=MPI.COMM_SELF,
        solver_run_scheme = "mpi_spawn_wrapper"
    )
    # output_base_dir = "/Users/k-yoshimi/Dropbox/PycharmProjects/abICS/test/openmx"
    output_base_dir = os.path.abspath(os.path.dirname(__file__))
    print(output_base_dir)
    energy_calculator.submit(structure=phys.structure, output_dir=os.path.join(output_base_dir, "test2"))
    # Check the coordinates of test1/met.dat# and test2/met.dat# are same.
