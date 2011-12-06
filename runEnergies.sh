#!/bin/sh

set -e

numParticles=300
density=0.01
LJcutoff=2.5
samples=300
measWait=0
measInterval=0.01
timestep=0.0005

./main $numParticles $density -m E -r -l $LJcutoff -s $samples -w $measWait -i $measInterval -t $timestep

octave --persist -q plotEnergies.m
