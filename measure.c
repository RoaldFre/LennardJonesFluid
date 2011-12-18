#include <stdlib.h>
#include "measure.h"
#include "system.h"

static void   calculatePairCorrelation(int bins, double *buf);
static void   addToDistanceBin(Particle *p1, Particle *p2, void *data);
static double pressureFromPairCorrelation(int bins, double *buf, double Temperature);

static void accumulatePairCorrelation(void);
static void accumulatePressure(void);
static void accumulateTemperature(void);

static void dumpAccumulatedPairCorrelation(FILE *stream);
static void dumpAccumulatedPressure(FILE *stream);
static void dumpEnergies(FILE *stream);
static void dumpPairCorrelation(FILE *stream, int bins, double *buf);

static void checkPairCorrelation(int bins, double *buf);


/* Buffer where we accumulate the pair correlation funcion measurements */
static double *accumPairCorrelationBuffer = NULL;

/* Bins of said buffer */
static int accumPairCorrelationBins = 100;

/* Amount of measurement samples made. */
static long samples = 0;


/* Instantaneous pressure */
double pressure(int bins)
{
	double *buf = calloc(bins, sizeof(double));
	calculatePairCorrelation(bins, buf);
	double P = pressureFromPairCorrelation(bins, buf, temperature());
	free(buf);
	return P;
}

static double pressureFromPairCorrelation(int bins, double *buf, double Temperature)
{
	/* P = rho kB T  -  1/6 rho^2 integral(4*pi*r^2 * r dphi/dr g(r) dr)
	 * where
	 * dphi/dr = -24 (2/r^13 - 1/r^7)
	 * combined, the integral becomes:
	 * integral(4*pi * (-24)*(2/r^10 - 1/r^4) * g(r) * dr) */
	double ws = config.numBox * config.boxSize;
	double rho = config.numParticles / (ws * ws * ws);

	double integral = 0;
	double dr = config.truncateLJ / bins;
	for (int i = 0; i < bins; i++) {
		double r = config.truncateLJ * (i + 1)/bins;
		double r2  = r*r;
		double r4  = r2*r2;
		double r8  = r4*r4;
		double r10  = r8*r2;
		
		integral += (2/r10 - 1/r4) * buf[i];
	}
	/* scale with the correct factors */
	integral *= 4 * M_PI * (-24) * dr;

	return rho * (Temperature - rho/6 * integral);
}

/* Check for consistency
 * NOTE: This only makes sense if the LJ-truncation length is sufficiently 
 * large (at least half the worldsize for periodic boundary conditions). 
 * This means everything reduces to O(n^2) because there will be no 
 * partitioning (one big box [or 2x2 boxes]) */
static void checkPairCorrelation(int bins, double *buf)
{
	/* We know:
	 * integral(rho * g(r) * 4*pi*r^2 dr) = N - 1
	 * So check if we recover this. */

	double ws = config.numBox * config.boxSize;
	double rho = config.numParticles / (ws * ws * ws);

	double integral = 0;
	double dr = config.truncateLJ / bins;
	for (int i = 0; i < bins; i++) {
		double r  = config.truncateLJ * (i + 1.0)/bins;
		integral += r*r * buf[i];
	}
	/* scale with the correct factors */
	integral *= rho * 4 * M_PI * dr;

	printf("Pair correlation integral: %f (num particles: %d)\n",
			integral, config.numParticles);
}


struct distanceBinData {
	int bins;
	double *buf;
};

/* Calculate the (discrete) pair correlation between r=0 and 
 * r=config.truncateLJ for the given number of bins in the given buffer. */
static void calculatePairCorrelation(int bins, double *buf)
{
	struct distanceBinData data;
	data.bins = bins;
	data.buf = buf;

	forEveryPair(&addToDistanceBin, &data);

	/* Rescale
	 * Don't forget: we only counted distinct pairs above, so we should 
	 * double everything to get 'every pair'.
	 * Then we should divide by N because we averaged over every pair, 
	 * while we actually should have considered one single particle 
	 * with every other particle (instead of *every* particle with 
	 * every other particle). */
	double ws = config.numBox * config.boxSize;
	double rho = config.numParticles / (ws * ws * ws);
	double dr = config.truncateLJ / bins;
	double factor = 2.0 / (rho * 4 * M_PI * dr * config.numParticles);
	for (int i = 0; i < bins; i++) {
		double r = config.truncateLJ * ((double) (i + 1.0) / (double) bins);
		buf[i] *= factor / (r*r);
	}
}

static void addToDistanceBin(Particle *p1, Particle *p2, void *data)
{
	struct distanceBinData *dbd = (struct distanceBinData *) data;
	int bins = dbd->bins;
	double *buf = dbd->buf;

	double Rmax = config.truncateLJ;
	Vec3 drVec = nearestImageVector(&p1->pos, &p2->pos);
	double dr = length(&drVec);
	if (dr >= Rmax)
		return;

	int i = bins * dr / Rmax;
	buf[i]++;
}

static void dumpPairCorrelation(FILE *stream, int bins, double *buf)
{
	for (int i = 0; i < bins; i++) {
		double r = config.truncateLJ * (i + 1)/bins;
		fprintf(stream, "%f\t%f\n", r, buf[i]);
	}
}

double temperature()
{
	return 2.0/3.0 * kineticEnergy() / config.numParticles;
}

static double lennardJonesPotential(Particle *p1, Particle *p2)
{
	assert(p1 != p2);
	assert(p1 != NULL  &&  p2 != NULL);

	Vec3 drVec = nearestImageVector(&p1->pos, &p2->pos);
	double dr = length(&drVec);

	assert(dr != 0);

	if (dr > config.truncateLJ)
		return 0;

	double dr3 = dr * dr * dr;	
	double dr6 = dr3 * dr3;	
	double dr12 = dr6 * dr6;	

	double dRc = config.truncateLJ;
	double dRc3 = dRc * dRc * dRc;
	double dRc6 = dRc3 * dRc3;
	double dRc12 = dRc6 * dRc6;

	return 4*(1/dr12 - 1/dr6 - (1/dRc12 - 1/dRc6));
}

static void potentialEnergyWorker(Particle *p1, Particle *p2, void *data)
{
	double *V = (double*) data;
	*V += lennardJonesPotential(p1, p2);
}

double potentialEnergy()
{
	double V = 0;
	forEveryPair(&potentialEnergyWorker, &V);
	return V;
}

double kineticEnergy()
{
	double totKE = 0;
	for (int i = 0; i < config.numParticles; i++)
	{
		double vsq;
		Particle *p = &world.parts[i];
		vsq = dot(&p->vel, &p->vel);
		totKE += vsq;
	}
	return totKE/2;
}

Vec3 momentum()
{
	Vec3 totP = {0, 0, 0};
	for (int i = 0; i < config.numParticles; i++)
	{
		Particle *p = &world.parts[i];
		add(&totP, &p->vel, &totP);
	}
	return totP;
}



/* GENERIC INTERFACE FOR MEASUREMENTS BY SAMPLING */

/* Will be called everytime we want to sample the system. Either dump the 
 * data immediately, or accumulate it to dump at the end.
 * Returns false if something went wrong. */
bool sampleMeasurement(FILE *stream)
{
	switch (config.measurement) {
	case NOTHING:
		break;
	case PRESSURE:
		accumulatePressure();
		break;
	case PAIR_CORRELATION:
		accumulatePairCorrelation();
		break;
	case ENERGIES:
		dumpEnergies(stream);
		break;
	default:
		fprintf(stderr, "Oops! Unknown measurement type! (type %d)\n", 
							config.measurement);
		break;
	}

	samples++;
	return true;
}

/* Will be called at the end of the measurement, to dump potentially 
 * accumulated data
 * Returns false if something went wrong. */
bool dumpMeasurement(FILE *stream)
{
	if (samples <= 0)
		return true; /* or maybe false, ... whatevs */

	switch (config.measurement) {
	case NOTHING:
	case ENERGIES:
		break;
	case PAIR_CORRELATION:
		dumpAccumulatedPairCorrelation(stream);
		break;
	case PRESSURE:
		dumpAccumulatedPressure(stream);
		break;
	default:
		fprintf(stderr, "Oops! Unknown measurement type! (type %d)\n", 
							config.measurement);
		return false;
	}
	return true;
}


/* INDIVIDUAL MEASUREMENTS */

static void accumulatePairCorrelation(void)
{
	int bins = accumPairCorrelationBins;
	if (accumPairCorrelationBuffer == NULL)
		accumPairCorrelationBuffer = calloc(bins, sizeof(double));

	double *buf = calloc(bins, sizeof(double));
	calculatePairCorrelation(bins, buf);
	for (int i = 0; i < bins; i++)
		accumPairCorrelationBuffer[i] += buf[i];
	free(buf);
}

static void dumpAccumulatedPairCorrelation(FILE *stream)
{
	if (accumPairCorrelationBuffer == NULL) {
		//fprintf(stderr, "accumPairCorrelationBuffer is NULL!");
		return;
	}
	int bins = accumPairCorrelationBins;
	for (int i = 0; i < bins; i++)
		accumPairCorrelationBuffer[i] /= samples;
	dumpPairCorrelation(stream, bins, accumPairCorrelationBuffer);
}

static void dumpEnergies(FILE *stream)
{
	double K = kineticEnergy();
	double V = potentialEnergy();

	sanityCheck();
	fprintf(stream, "%13f\t%13f\t%13f\t%13f\n", sim_time, K+V, K, V);
}

/* Accumulate temperature measurements. */
static double accumTemperature = 0;

static void accumulateTemperature(void)
{
	accumTemperature += temperature();
}

static void accumulatePressure(void)
{
	accumulatePairCorrelation();
	accumulateTemperature();
}

static void dumpAccumulatedPressure(FILE *stream)
{
	for (int i = 0; i < accumPairCorrelationBins; i++)
		accumPairCorrelationBuffer[i] /= samples;

	/*
	checkPairCorrelation(accumPairCorrelationBins,
			accumPairCorrelationBuffer);
	*/

	fprintf(stream, "%f", pressureFromPairCorrelation(
				accumPairCorrelationBins,
				accumPairCorrelationBuffer,
				accumTemperature / samples));
}

/* Check if the physics still make sense. */
bool physicsCheck()
{
	Vec3 P = momentum();
	double Plength = length(&P);
	if (Plength > 0.0001) {
		fprintf(stderr, "\nMOMENTUM CONSERVATION VIOLATED! "
				"Total momentum: |P| = %f\n", Plength);
		return false;
	}

	return true;
}
