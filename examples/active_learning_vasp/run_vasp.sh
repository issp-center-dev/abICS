# specify program
prog="srun --exclusive --mem-per-cpu=1840 -n 32 -c 2 -N 1 vasp_gam" 

cd $1
$prog > stdout
