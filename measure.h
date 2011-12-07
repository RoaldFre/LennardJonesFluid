#ifndef _MEASURE_H_
#define _MEASURE_H_

#include <stdio.h>
#include <stdbool.h>
#include "vmath.h"

enum measurementType {
	NOTHING,
	PAIR_CORRELATION,
	ENERGIES,
	PRESSURE
};

Vec3   momentum(void);
double temperature(void);
double kineticEnergy(void);
double potentialEnergy(void);
double pressure(int bins);

bool sampleMeasurement(FILE *stream);
bool dumpMeasurement(FILE *stream);
bool physicsCheck(void);

#endif 
