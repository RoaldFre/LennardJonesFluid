#!/bin/sh

#WARNING: results are invalid for radi larger than half(?) the boxsize,
#-> adjust the lennard jones cut off accordingly!

set -e

./main 100 0.01 -r -l 3 -s 300 -w 20 -i 0.2 -t 0.0005

octave --persist plotData.m 
