#!/bin/sh

#WARNING: results are invalid for radi larger than half(?) the boxsize,
#-> adjust the lennard jones cut off accordingly!

set -e

numParticles=200
density=0.8
LJcutoff=-1
samples=100
measWait=10
measInterval=0.2
timestep=0.005

./main $numParticles $density -m c -r -l $LJcutoff -s $samples -w $measWait -i $measInterval -t $timestep

octave --persist -q plotPairCorrelation.m
