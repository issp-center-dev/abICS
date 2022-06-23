mkdir -p pseudo
for file in \
Al.pbe-nl-kjpaw_psl.1.0.0.UPF \
Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF \
O.pbe-n-kjpaw_psl.1.0.0.UPF
do
wget https://pseudopotentials.quantum-espresso.org/upf_files/$file
mv $file pseudo
done
