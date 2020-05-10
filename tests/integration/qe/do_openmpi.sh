OMP_NUM_THREADS=1 mpiexec -np 2 --oversubscribe --mca orte_abort_on_non_zero_status 0 python ./fulltest.py
