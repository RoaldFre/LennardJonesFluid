#!/bin/sh

set -e

numParticles=300
LJcutoff=2.5
samples=200
measWait=0
measInterval=0.05
timestep=0.01
coupling=1

temperature=2
initialTemp="0.001"

for density in `seq 0.2 0.2 0.8`
do
	./main $numParticles $density -m E -r -l $LJcutoff -s $samples -w $measWait -i $measInterval -t $timestep -T $temperature -c $coupling
	densityNoPoint=`echo $density | tr . p`
	filename="energies-$densityNoPoint"
	info="\$N\$ = $numParticles, \$\\\\rho\$ = $density, \$\Delta t\$ = $timestep, \$\\\\tau\$ = $coupling, \$T\$ = $temperature, \$T_0 = T\$, \$\\\\tau_\\\\mathrm{sample}\$ = $measInterval, \$N_\\\\mathrm{samples}\$ = $samples"

	caption="Potential, kinetic and total energy, starting from equilibrium temperature. Parameters: $info."
	octave -q --eval "plotEnergies('$filename', '$caption')"

	./main $numParticles $density -m E -r -l $LJcutoff -s $samples -w $measWait -i $measInterval -t $timestep -T $temperature -I $initialTemp -c $coupling
	densityNoPoint=`echo $density | tr . p`
	filename="energies-T0-$densityNoPoint"
	info="\$N\$ = $numParticles, \$\\\\rho\$ = $density, \$\Delta t\$ = $timestep, \$\\\\tau\$ = $coupling, \$T\$ = $temperature, \$T_0\$ = $initialTemp, \$\\\\tau_\\\\mathrm{sample}\$ = $measInterval, \$N_\\\\mathrm{samples}\$ = $samples"
	caption="Relaxation of potential, kinetic and total energy starting from low temperature. Parameters: $info."
	octave -q --eval "plotEnergies('$filename', '$caption')"
done

