#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#define MEM_ERROR(n) (fprintf(stderr, "Not enough memory for %ld bytes\n", (long) n))

#include <stdbool.h>
#include "vmath.h"

typedef struct particle
{
	Vec3 pos; /* Position */
	Vec3 vel; /* Velocity */
	Vec3 acc; /* Acceleration */
	struct particle *prev, *next; /* Previous/Next particle in box */
} Particle;

typedef struct box
{
	Particle *p;
	int n;
} Box;

typedef struct world
{
	Box *grid;
	Particle *parts;
} World;

typedef struct config
{
	bool render;
	int verbose;
	double boxSize; /* The same for every dimension */
	int numBox;
	int numParticles;
	long   measureSamples;  /* Total number of samples to accumulate */
	double measureInterval; /* Time between samples */
	double measureWait;     /* Time to wait before starting measurement */
	long   renderSteps;     /* Physics steps between rendering frames. */
	double timeStep;
	double truncateLJ; /* Radius at which L-J potential gets truncated. */
	double temperature;    /* Desired temperature. */
	double thermostatTemp; /* Thermostat temperature. */
	double thermostatTau;  /* Thermostat relaxation time. */
	double radius; /* The radius of the particles to render */
} Config;

int main(int argc, char ** argv);
bool allocWorld(void);
void fillWorld(void);
void freeWorld(void);
void dumpWorld(void);
void stepWorld(void);
bool sanityCheck(void);
bool physicsCheck(void);
void dumpStats(void);
void dumpPairCorrelation(FILE *stream, int bins, double *buf);
void dumpPairCorrelationm(FILE *stream, int bins, double *buf);
void accumulatePairCorrelation(void);
void dumpAccumulatedPairCorrelation(FILE *stream);
void sampleMeasurement(void);
void dumpMeasurement(FILE *stream);

extern World world;
extern Config config;

#endif
