#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#define MEM_ERROR(n) (fprintf(stderr, "Not enough memory for %ld bytes\n", (long) n))

#include <stdbool.h>
#include "vmath.h"
#include "measure.h"

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
	enum measurementType measurement;
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
	double truncateLJ;      /* Radius at which L-J potential gets truncated. */
	double initialTemp;     /* Initial temperature. */
	double thermostatTemp;  /* Thermostat temperature. */
	double thermostatTau;   /* Thermostat relaxation time. */
	double radius; /* The radius of the particles to render */
} Config;

bool allocWorld(void);
void fillWorld(void);
void freeWorld(void);
void stepWorld(void);
void dumpWorld(void);
void dumpStats(void);
bool sanityCheck(void);
void forEveryPair(void (*f)(Particle *p1, Particle *p2, void *data), void *data);
Vec3 nearestImageVector(Vec3 *v1, Vec3 *v2); 

extern World world;
extern Config config;

/* Current time in the simulation */
extern double sim_time;

#endif
