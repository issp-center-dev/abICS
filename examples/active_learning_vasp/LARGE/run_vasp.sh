# specify program
prog="srun --exclusive --mem-per-cpu=1840 -n 128 -c 2 -N 2 /home/k0306/k030600/src/vasp.5.4.4.pl2/bin/vasp_gam" 

cd $1
$prog > stdout
