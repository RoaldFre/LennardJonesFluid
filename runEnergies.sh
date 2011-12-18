#!/bin/sh

set -e

numParticles=300
LJcutoff=2.5
samples=200
measWait=0
measInterval=0.01
timestep=0.01

temperature=2
initialTemp="1e-5"

for density in `seq 0.2 0.2 0.8`
do
	./main $numParticles $density -m E -r -l $LJcutoff -s $samples -w $measWait -i $measInterval -t $timestep -T $temperature
	densityNoPoint=`echo $density | tr . p`
	filename="energies-$densityNoPoint"
	info="\$N\$ = $numParticles, \$\\\\rho\$ = $density, \$\Delta t\$ = $timestep, \$\\\\tau\$ = \$10^3 \Delta t\$, \$T\$ = $temperature"
	caption="Potential, kinetic and total energy. Parameters: $info."
	octave -q --eval "plotEnergies('$filename', '$caption')"

	./main $numParticles $density -m E -r -l $LJcutoff -s $samples -w $measWait -i $measInterval -t $timestep -T $temperature -I $initialTemp
	densityNoPoint=`echo $density | tr . p`
	filename="energies-T0-$densityNoPoint"
	info="\$N\$ = $numParticles, \$\\\\rho\$ = $density, \$\Delta t\$ = $timestep, \$\\\\tau\$ = \$10^3 \Delta t\$, \$T\$ = $temperature, \$T_0\$ = $initialTemp"
	caption="Potential, kinetic and total energy. Parameters: $info."
	octave -q --eval "plotEnergies('$filename', '$caption')"
done

