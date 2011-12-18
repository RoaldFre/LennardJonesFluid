#!/bin/sh

set -e

numParticles=300
LJcutoff=2.5
samples=50
measWait=10
measInterval=0.1
timestep=0.005
temperature=2

file=/tmp/pressure.txt

echo -n "" > $file

for density in `seq 0.01 0.03 1`
do
	./main $numParticles $density -m P -r -l $LJcutoff -s $samples -w $measWait -i $measInterval -t $timestep -T $temperature

	echo -n -e "$density\t" >> $file
	cat /tmp/data.txt >> $file
	echo "" >> $file
done



info="\$N\$ = $numParticles, \$\Delta t\$ = $timestep, \$\\\\tau\$ = \$10^3 \Delta t\$, \$T\$ = $temperature, \$\\\\tau_\\\\mathrm{relax}\$ = $measWait, \$\\\\tau_\\\\mathrm{sample}\$ = $measInterval, \$N_\\\\mathrm{samples}\$ = $samples"
caption="Pressure in function of density. Parameters: $info."
filename="pressure"
octave -q --eval "plotPressure('$filename', '$caption', $temperature)"
