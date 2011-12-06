#!/bin/sh

set -e

numParticles=200
LJcutoff=-1
samples=50
measWait=10
measInterval=0.15
timestep=0.005

file=/tmp/pressure.txt

echo -n "" > $file

for density in `seq 0.01 0.05 1`
do
	./main $numParticles $density -m P -r -l $LJcutoff -s $samples -w $measWait -i $measInterval -t $timestep

	echo -n -e "$density\t" >> $file
	cat /tmp/data.txt >> $file
	echo "" >> $file
done
	




octave --persist -q plotPressure.m
