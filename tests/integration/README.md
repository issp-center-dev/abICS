# Integration test

These take ten or more minutes to finish.

## QE and QE-restart
1. Put binary `pw.x` into `~/opt/qe/bin`
2. Get the following pseudopotentials
```
Al.pbe-nl-kjpaw_psl.1.0.0.UPF
Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF
O.pbe-n-kjpaw_psl.1.0.0.UPF
```
and put them into `~/opt/qe/pot`
3. Invoke `do_openmpi.sh` or `do_mpich.sh`

## OpenMX
1. Put binary `openmx` into `~/opt/openmx/bin`
2. Put datafiles `DFT_DATA19` into `~/opt/openmx/`
3. Invoke `do_openmpi.sh` or `do_mpich.sh`
