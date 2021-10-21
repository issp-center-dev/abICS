# specify program

source /home/issp/materiapps/intel/espresso/espressovars.sh

pwd
srun --exclusive --mem-per-cpu=1840 -n 32 -c 1 -N 1 pw.x -in $1/scf.in
echo "$1" finished
