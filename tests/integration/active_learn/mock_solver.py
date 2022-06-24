from abics.applications.latgas_abinitio_interface.mocksolver import MockSolver

def main(output_dir):
    solver = MockSolver()
    solver.calc_energy("", output_dir)

if __name__ == "__main__":
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    main(output_dir)
