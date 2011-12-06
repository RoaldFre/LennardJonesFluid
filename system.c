#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "vmath.h"
#include "system.h"

World world;
Config config;

static void lennardJonesForce(Particle *p1, Particle *p2);
static void collideWalls(int ix, int iy, int iz);
static void reboxParticles(void);
static Box *boxFromParticle(const Particle *p);
static Box *boxFromIndex(int nx, int ny, int nz);
static Box *boxFromNonPeriodicIndex(int nx, int ny, int nz);
static void addToBox(Particle *p, Box *b);
static void removeFromBox(Particle *p, Box *from);
static void calculateAcceleration(void);
static void calculateAccelerationBruteForce(void);
static void verlet(void);
static double randNorm(void);
static Vec3 momentum(void);
static double kineticEnergy(void);
static double temperature(void);
static void thermostat(void);
static double potentialEnergy(void);
static double potentialEnergyBruteForce(void);
static void calculatePairCorrelation(int bins, double *buf);
static void addPairCorrelation(Particle *p1, Particle *p2, int bins, double *buf);
static Vec3 nearestImageVector(Vec3 *v1, Vec3 *v2); 

/* UNDER THE BONNET STUFF */

/* Put particles back in their correct boxes. */
static void reboxParticles(void)
{
	int nb = config.numBox;

	for (int i = 0; i < nb*nb*nb; i++) {
		Box *currentBox = &world.grid[i];
		Particle *p = currentBox->p;

		if (p == NULL)
			continue; /* Empty box */

		int n = currentBox->n; /* Save original number, it can 
					  decrease during the loop if we 
					  swap out particles! */
		for (int j = 0; j < n; j++) {
			/* Since p might be removed from the current box, 
			 * we keep a pointer to its successor. */
			Particle *next = p->next;

			Box *correctBox = boxFromParticle(p);
			if (currentBox != correctBox) {
				removeFromBox(p, currentBox);
				addToBox(p, correctBox);
			}

			p = next;
		}
	}
	assert(sanityCheck());
}

static Box *boxFromParticle(const Particle *p)
{
	int nx, ny, nz;
	double bs = config.boxSize;

	assert(p != NULL);
	assert(!isnan(p->pos.x) && !isnan(p->pos.y) && !isnan(p->pos.z));
	assert(0 <= p->pos.x  &&  p->pos.x < bs * config.numBox);
	assert(0 <= p->pos.y  &&  p->pos.y < bs * config.numBox);
	assert(0 <= p->pos.z  &&  p->pos.z < bs * config.numBox);

	nx = p->pos.x / bs;
	ny = p->pos.y / bs;
	nz = p->pos.z / bs;

	return boxFromIndex(nx, ny, nz);
}

static Box *boxFromNonPeriodicIndex(int ix, int iy, int iz)
{
	int nb = config.numBox;

	ix = ix % nb;
	if (ix < 0) ix += nb;

	iy = iy % nb;
	if (iy < 0) iy += nb;

	iz = iz % nb;
	if (iz < 0) iz += nb;

	return boxFromIndex (ix, iy, iz);
}

static Box *boxFromIndex(int ix, int iy, int iz)
{
	int nb = config.numBox;
	
	assert(0 <= ix && ix < nb);
	assert(0 <= iy && iy < nb);
	assert(0 <= iz && iz < nb);

	return world.grid + ix*nb*nb + iy*nb + iz;
}

static void removeFromBox(Particle *p, Box *b)
{
	assert(p != NULL);
	assert(b->n != 0);

	if (b->n == 1) {
		assert(p->prev == p);
		assert(p->next == p);
		b->p = NULL;
	} else {
		assert(p->prev->next == p);
		assert(p->next->prev == p);
		p->prev->next = p->next;
		p->next->prev = p->prev;

		if (b->p == p)
			b->p = p->next;
	}

	p->prev = NULL;
	p->next = NULL;

	b->n--;
}

static void addToBox(Particle *p, Box *b)
{
	assert(p->prev == NULL);
	assert(p->next == NULL);

	if (b->p == NULL) {
		assert(b->n == 0);
		b->p = p;
		p->prev = p;
		p->next = p;
	} else {
		assert(b->n > 0);
		p->next = b->p;
		p->prev = b->p->prev;
		p->prev->next = p;
		p->next->prev = p;
	}

	b->n++;
}

/* Precondition: config MUST be valid
 * Allocates */
bool allocWorld()
{
	int nb = config.numBox;

	world.parts = calloc(config.numParticles, sizeof(Particle));
	if (world.parts == NULL) {
		MEM_ERROR(config.numParticles * sizeof(Particle));
		return false;
	}

	world.grid = calloc(nb * nb * nb, sizeof(*(world.grid)));
	if (world.grid == NULL) {
		MEM_ERROR(nb * nb * nb * sizeof(*world.grid));
		return false;
	}
	return true;
}

/* Place particles at cubic lattice in the center, a distance 
 * config.truncateLJ apart.
 * If that is not possible (because the world is too small), the particles 
 * are placed on a uniform lattice filling the entire world.
 *
 * Initial velocities are sampled from a normal distribution with variance 
 * 'config.temperature'. Velocities are then shifted to set the total 
 * momentum to zero. */
void fillWorld()
{
	Particle *ps = world.parts;
	double worldSize = config.numBox * config.boxSize;
	double stdev = sqrt(config.temperature); // TODO factor 1/3 (3D) somethere?
	int nPerDim = ceil(cbrt(config.numParticles));

	double spacing = config.truncateLJ;

	/* If it's too close, pack evenly as best as we can */
	//if (spacing * nPerDim > worldSize) // always do this, makes much more sense
		spacing = worldSize / nPerDim;

	Vec3 totVel = {0, 0, 0}, avgVel;

	for (int i = 0; i < config.numParticles; i++) {
		int ix = i % nPerDim;
		int iy = (i / nPerDim) % nPerDim;
		int iz = i / (nPerDim * nPerDim);

		ps[i].pos.x = ix*spacing + (worldSize - spacing*(nPerDim-1))/2;
		ps[i].pos.y = iy*spacing + (worldSize - spacing*(nPerDim-1))/2;
		ps[i].pos.z = iz*spacing + (worldSize - spacing*(nPerDim-1))/2;

		// TODO: check before to make sure we don't need this, but 
		//here for now just to be safe
		periodic(worldSize, &ps[i].pos, &ps[i].pos);

		ps[i].vel.x = stdev * randNorm();
		ps[i].vel.y = stdev * randNorm();
		ps[i].vel.z = stdev * randNorm();
		add(&totVel, &ps[i].vel, &totVel);
		
		ps[i].acc.x = 0; ps[i].acc.y = 0; ps[i].acc.z = 0;

		Box *box = boxFromParticle(&ps[i]);
		addToBox(&ps[i], box);
	}
	scale(&totVel, 1.0 / config.numParticles, &avgVel);
	for (int i = 0; i < config.numParticles; i++)
		sub(&ps[i].vel, &avgVel, &ps[i].vel);

	assert(sanityCheck());
	printf("Initial situation:\n");
	dumpStats();
	printf("\n");
}

/* Returns a number sampled from a standard normal distribution. */
static double randNorm()
{
	/* Box-Muller transform */
	double u1 = ((double) rand()) / RAND_MAX;
	double u2 = ((double) rand()) / RAND_MAX;

	return sqrt(-2*log(u1)) * cos(2*M_PI*u2);
}

void freeWorld()
{
	free(world.grid);
	free(world.parts);

	return;
}

void dumpWorld()
{
	for (int i = 0; i < config.numParticles; i++) {
		const Particle *p = &world.parts[i];
		printf("%d\t", i);
		printVector(&p->pos);
		printVector(&p->vel);
		printf("\n");
	}

	return;
}



/* PHYSICS */

static void verlet()
{
	double dt = config.timeStep;
	double worldSize = config.numBox * config.boxSize;

	/* Velocity Verlet */
	for (int i = 0; i < config.numParticles; i++) {
		Particle *p = &world.parts[i];
		Vec3 tmp;

		/* vel(t + dt/2) = vel(t) + acc(t)*dt/2 */
		scale(&p->acc, dt/2, &tmp);
		add(&p->vel, &tmp, &p->vel);

		assert(!isnan(p->vel.x) && !isnan(p->vel.y) && !isnan(p->vel.z));

		/* pos(t + dt) = pos(t) + vel(t + dt/2)*dt */
		scale(&p->vel, dt, &tmp);
		add(&p->pos, &tmp, &p->pos);

		/* Periodic boundary conditions */
		periodic(worldSize, &p->pos, &p->pos);
	}
	reboxParticles(); /* Move particles to correct box if they escaped */
	calculateAcceleration(); /* acc(t + dt) */
	//calculateAccelerationBruteForce(); /* acc(t + dt) */
	for (int i = 0; i < config.numParticles; i++) {
		Particle *p = &world.parts[i];
		Vec3 tmp;

		/* vel(t + dt) = vel(t + dt/2) + acc(t + dt)*dt/2 */
		scale(&p->acc, dt/2, &tmp);
		add(&p->vel, &tmp, &p->vel);
	}
}

static void calculateAcceleration()
{
	int i, j, ix, iy, iz, dix, diy, diz, n1, n2;
	int nb = config.numBox;

	assert(sanityCheck());

	if (nb < 3) {
		calculateAccelerationBruteForce();
		/* Reason: see comment below */
		return;
	}

	/* Reset acceleration */
	for (i = 0; i < config.numParticles; i++) {
		world.parts[i].acc.x = 0;
		world.parts[i].acc.y = 0;
		world.parts[i].acc.z = 0;
	}
	
	/* Calculate acceleration */
	for (ix = 0; ix < nb; ix++)
	for (iy = 0; iy < nb; iy++)
	for (iz = 0; iz < nb; iz++) {
		/* Calculate acceleration of all particles in this box */
		Box *box = boxFromIndex(ix, iy, iz);
		Particle *p = box->p;

		/* Calculate acceleration of the i'th particles 'p' in the 
		 * box 'box' */
		n1 = box->n;
		for (i = 0; i < n1; i++) {
			/* Loop over every pair of particles in the box of p */
			Particle *p2 = p->next;
			for (j = i + 1; j < n1; j++) {
				lennardJonesForce(p, p2);
				p2 = p2->next;
			}
			/* Loop over particles in adjacent boxes to the box 
			 * of p. We need a total ordering on the boxes so 
			 * we don't check the same box twice. We use the 
			 * pointer value for this.
			 * However, due to periodic boundary conditions, 
			 * this ONLY works when there are AT LEAST 3 boxes 
			 * in each dimension! */
			for (dix = -1; dix <= 1; dix++)
			for (diy = -1; diy <= 1; diy++)
			for (diz = -1; diz <= 1; diz++) {
				Box *b = boxFromNonPeriodicIndex(
						ix+dix, iy+diy, iz+diz);
				if (b <= box)
					continue;
					/* if b == box: it's our own box!
					 * else: only check boxes that have 
					 * a strictly smaller pointer value 
					 * to avoid double work. */
				p2 = b->p;
				n2 = b->n;
				for (j = 0; j < n2; j++) {
					lennardJonesForce(p, p2);
					p2 = p2->next;
				}
			}

			p = p->next;
		}
	}
}

static void calculateAccelerationBruteForce()
{
	/* Reset acceleration */
	for (int i = 0; i < config.numParticles; i++) {
		world.parts[i].acc.x = 0;
		world.parts[i].acc.y = 0;
		world.parts[i].acc.z = 0;
	}
	
	for (int i = 0; i < config.numParticles; i++) {
		Particle *p1 = &world.parts[i];
		for (int j = i + 1; j < config.numParticles; j++) {
			Particle *p2 = &world.parts[j];
			lennardJonesForce(p1, p2);
		}
	}
}

/* Returns the shortest vector that points from v1 to v2, taking into 
 * account the periodic boundary conditions. The particle MUST be within 
 * the correct bounds. */
static Vec3 nearestImageVector(Vec3 *v1, Vec3 *v2)
{
	Vec3 diff;
	double L = config.numBox * config.boxSize;
	/* +5*L/2 instead of simply +L/2 to make sure that the first 
	 * argument of fmod is positive. Otherwise, this won't work! */
	diff.x = fmod(v2->x - v1->x + 5*L/2, L)  -  L/2;
	diff.y = fmod(v2->y - v1->y + 5*L/2, L)  -  L/2;
	diff.z = fmod(v2->z - v1->z + 5*L/2, L)  -  L/2;
	return diff;
}

/* Adds the acceleration due to the LJ force to the acceleration of the 
 * given particles. */
static void lennardJonesForce(Particle *p1, Particle *p2)
{
	assert(p1 != p2);
	assert(p1 != NULL  &&  p2 != NULL);

	Vec3 Fi;
	Vec3 drVec = nearestImageVector(&p1->pos, &p2->pos);
	double dr = length(&drVec);

	assert(dr != 0);

	if (dr > config.truncateLJ)
		return;

	double dr2 = dr*dr;	
	double dr3 = dr*dr*dr;	
	double dr6 = dr3*dr3;	
	double dr8 = dr6*dr2;	
	double dr12 = dr6*dr6;	
	double dr14 = dr12*dr2;	


	//assert(fabs(-24*(2/dr12 - 1/dr6)) < 1000000);

	scale(&drVec, -24*(2/dr14 - 1/dr8), &Fi);
	add(&p1->acc, &Fi, &p1->acc);
	sub(&p2->acc, &Fi, &p2->acc);
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

static void collideWalls(int ix, int iy, int iz)
{
	int i;
	int nb = config.numBox;
	Particle *p;
	Box *b = boxFromIndex(ix, iy, iz);
	double worldSize = config.boxSize * config.numBox;

	p = b->p;
	for (i = 0; i < b->n; i++)
	{
		if (ix == 0 && p->pos.x < 0)
		{
			p->vel.x = -p->vel.x;
			p->pos.x = -p->pos.x;
		} else if (ix == nb && p->pos.x > worldSize)
		{
			p->vel.x = -p->vel.x;
			p->pos.x = worldSize - (p->pos.x - worldSize);
		}

		if (iy == 0 && p->pos.y < 0)
		{
			p->vel.y = -p->vel.y;
			p->pos.y = -p->pos.y;
		} else if (iy == nb && p->pos.y > worldSize)
		{
			p->vel.y = -p->vel.y;
			p->pos.y = worldSize - (p->pos.y - worldSize);
		}

		if (iz == 0 && p->pos.z < 0)
		{
			p->vel.z = -p->vel.z;
			p->pos.z = -p->pos.z;
		} else if (iz == nb && p->pos.z > worldSize)
		{
			p->vel.z = -p->vel.z;
			p->pos.z = worldSize - (p->pos.z - worldSize);
		}
	
		p = p->next;
	}
}

static void thermostat()
{
	if (config.thermostatTau <= 0)
		return;

	/* Mass and Boltzmann constant are 1 */ 
	double Tk  = temperature();
	double T0  = config.temperature;
	double dt  = config.timeStep;
	double tau = config.thermostatTau;
	double lambda = sqrt(1 + dt/tau * (T0/Tk - 1));

	for (int i = 0; i < config.numParticles; i++) {
		Particle *p = &world.parts[i];
		scale(&p->vel, lambda, &p->vel);
	}
}

void stepWorld(void)
{
	assert(sanityCheck());
	verlet();
	assert(sanityCheck());
	thermostat();
	assert(sanityCheck());
}


/* PHYSICS MEASUREMENTS */

static double pressure(int bins)
{
	/* P = rho kB T  -  1/6 rho^2 integral(4*pi*r^3 dphi/dr g(r) dr)
	 * where
	 * dphi/dr = -24 (2/r^13 - 1/r^7)
	 * combined, the integral becomes:
	 * integral(4*pi * (-24)*(2/r^10 - 1/r^4) * g(r) * dr) */
	double T = temperature();
	double ws = config.numBox * config.boxSize;
	double rho = config.numParticles / (ws * ws * ws);

	double *buf = calloc(bins, sizeof(double));
	calculatePairCorrelation(bins, buf);

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

	free(buf);

	return rho * (T - rho/6 * integral);
}

/* Calculate the discrete pair correlation between r=0 and 
 * r=config.truncateLJ for the given number of bins in the given buffer. */
static void calculatePairCorrelation(int bins, double *buf)
{
	int nb = config.numBox;
	double ws = nb * config.boxSize;
	double rho = config.numParticles / (ws * ws * ws);

	assert(sanityCheck());

	for (int i = 0; i < bins; i++) {
		buf[i] = 0;
	}

	if (nb < 3) {
		/* Brute force it */
		for (int i = 0; i < config.numParticles; i++) {
			Particle *p1 = &world.parts[i];
			for (int j = i + 1; j < config.numParticles; j++) {
				Particle *p2 = &world.parts[j];
				addPairCorrelation(p1, p2, bins, buf);
			}
		}
	} else {
		/* Use the boxes */
		for (int ix = 0; ix < nb; ix++)
		for (int iy = 0; iy < nb; iy++)
		for (int iz = 0; iz < nb; iz++) {
			/* Loop over all boxes */
			Box *box = boxFromIndex(ix, iy, iz);
			Particle *p = box->p;

			for (int i = 0; i < box->n; i++) { /* i'th particle in box */
				Particle *p2 = p->next;
				for (int j = i + 1; j < box->n; j++) {
					addPairCorrelation(p, p2, bins, buf);
					p2 = p2->next;
				}
				for (int dix = -1; dix <= 1; dix++)
				for (int diy = -1; diy <= 1; diy++)
				for (int diz = -1; diz <= 1; diz++) {
					Box *b = boxFromNonPeriodicIndex(
							ix+dix, iy+diy, iz+diz);
					if (b <= box)
						continue;
					p2 = b->p;
					for (int j = 0; j < b->n; j++) {
						addPairCorrelation(p, p2, bins, buf);
						p2 = p2->next;
					}
					assert(p2 == b->p);
				}

				p = p->next;
			}

			assert(p == box->p);
		}
	}

	/* Rescale
	 * Don't forget: we only counted distinct pairs above, so we should 
	 * double everything to get 'every pair'. */
	double dr = config.truncateLJ / bins;
	double factor = 2.0 / (rho * 4 * M_PI * dr);
	for (int i = 0; i < bins; i++) {
		double r = config.truncateLJ * ((double) (i + 1) / (double) bins);
		buf[i] *= factor / (r*r);
	}
}

static void addPairCorrelation(Particle *p1, Particle *p2, int bins, double *buf)
{
	double Rmax = config.truncateLJ;
	Vec3 drVec = nearestImageVector(&p1->pos, &p2->pos);
	double dr = length(&drVec);
	if (dr >= Rmax)
		return;

	int i = bins * dr / Rmax;
	buf[i]++;
}

/* Dump pair correlation in GNU Octave / Matlab format */
void dumpPairCorrelation(FILE *stream, int bins, double *buf)
{
	for (int i = 0; i < bins; i++) {
		double r = config.truncateLJ * (i + 1)/bins;
		fprintf(stream, "%f\t%f\n", r, buf[i]);
	}
}

/* Dump pair correlation in GNU Octave / Matlab format */
void dumpPairCorrelationm(FILE *stream, int bins, double *buf)
{
	double r;
	r = config.truncateLJ * (1)/bins;
	fprintf(stream, "[%f, %f;\n", r, buf[0]);
	for (int i = 1; i < (bins - 1); i++) {
		r = config.truncateLJ * (i + 1)/bins;
		fprintf(stream, " %f, %f;\n", r, buf[i]);
	}
	r = config.truncateLJ;
	fprintf(stream, " %f, %f];\n", r, buf[bins - 1]);
}

static double temperature()
{
	return 2.0/3.0 * kineticEnergy() / config.numParticles;
}

static double potentialEnergy()
{
	double totV = 0;
	int nb = config.numBox;
	int n1, n2;

	assert(sanityCheck());

	if (nb < 3)
		return potentialEnergyBruteForce();
		/* Reason: see comment below */

	for (int ix = 0; ix < nb; ix++)
	for (int iy = 0; iy < nb; iy++)
	for (int iz = 0; iz < nb; iz++) {
		/* Loop over all boxes */
		Box *box = boxFromIndex(ix, iy, iz);
		Particle *p = box->p;

		n1 = box->n;
		for (int i = 0; i < n1; i++) { /* i'th particle in box */
			Particle *p2 = p->next;
			for (int j = i + 1; j < n1; j++) {
				totV += 2*lennardJonesPotential(p, p2);
				p2 = p2->next;
			}
			/* Loop over particles in adjacent boxes to the box 
			 * of p. We need a total ordering on the boxes so 
			 * we don't check the same box twice. We use the 
			 * pointer value for this.
			 * However, due to periodic boundary conditions, 
			 * this ONLY works when there are AT LEAST 3 boxes 
			 * in each dimension! */
			for (int dix = -1; dix <= 1; dix++)
			for (int diy = -1; diy <= 1; diy++)
			for (int diz = -1; diz <= 1; diz++) {
				Box *b = boxFromNonPeriodicIndex(
						ix+dix, iy+diy, iz+diz);
				if (b == box)
					continue;
					/* if b == box: it's our own box!
					 * else: only check boxes that have 
					 * a strictly smaller pointer value 
					 * to avoid double work. */
				p2 = b->p;
				n2 = b->n;
				for (int j = 0; j < n2; j++) {
					totV += lennardJonesPotential(p, p2);
					p2 = p2->next;
				}
				assert(p2 == b->p);
			}

			p = p->next;
		}

		assert(p == box->p);
	}

	return totV/2;
}

static double potentialEnergyBruteForce()
{
	double totV = 0;

	for (int i = 0; i < config.numParticles; i++) {
		Particle *p1 = &world.parts[i];
		for (int j = i + 1; j < config.numParticles; j++) {
			Particle *p2 = &world.parts[j];
			totV += lennardJonesPotential(p1, p2);
		}
	}

	return totV;
}

static double kineticEnergy()
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

static Vec3 momentum()
{
	Vec3 totP = {0, 0, 0};
	for (int i = 0; i < config.numParticles; i++)
	{
		Particle *p = &world.parts[i];
		add(&totP, &p->vel, &totP);
	}
	return totP;
}

void dumpStats()
{
	double K = kineticEnergy();
	double V = potentialEnergy();
	//double T = temperature();
	double T = 2.0/3.0 * kineticEnergy() / config.numParticles;
	double P = pressure(200);	

	sanityCheck();
	printf("E = %13f \tK = %13f \tV = %13f \tT = %13f \tP = %13f\n", K+V, K, V, T, P);
	//printf("dV = %13f\n", V - potentialEnergyBruteForce());
}



/* GENERIC INTERFACE FOR MEASUREMENTS BY SAMPLING */

static long samples = 0;

/* Will be called everytime we want to sample the system */
void sampleMeasurement()
{
	accumulatePairCorrelation();
	samples++;
}

/* Will be called at the end of the measurement, to dump the data */
void dumpMeasurement(FILE *stream)
{
	dumpAccumulatedPairCorrelation(stream);
}


/* INDIVIDUAL MEASUREMENTS BY SAMPLING */

static double *pairCorrelationBuffer;
static int pairCorrelationBins = 300;

void accumulatePairCorrelation()
{
	int bins = pairCorrelationBins;
	if (pairCorrelationBuffer == NULL)
		pairCorrelationBuffer = calloc(bins, sizeof(double));

	double *buf = calloc(bins, sizeof(double));
	calculatePairCorrelation(bins, buf);
	for (int i = 0; i < bins; i++)
		pairCorrelationBuffer[i] += buf[i];
	free(buf);
}

void dumpAccumulatedPairCorrelation(FILE *stream)
{
	if (pairCorrelationBuffer == NULL) {
		fprintf(stderr, "pairCorrelationBuffer is NULL!");
		return;
	}
	int bins = pairCorrelationBins;
	for (int i = 0; i < bins; i++)
		pairCorrelationBuffer[i] /= samples;
	dumpPairCorrelation(stream, bins, pairCorrelationBuffer);
}



/* CHECKS */

/* Check whether internal structure is still consistent. */
bool sanityCheck()
{
	int i, j, nParts1, nParts2;
	const Particle *p1;
	nParts1 = 0;
	nParts2 = 0;
	int nb = config.numBox;
	bool OK = true;

	/* Make sure particle are not too close together and that their linked
	 * list is consistent */
	for (i = 0; i < config.numParticles; i++)
	{
		p1 = &world.parts[i];
		if (p1->next->prev != p1 || p1->prev->next != p1) {
			fprintf(stderr, "%p is in a borked list\n",
					(const void *) p1);
			OK = false;
		}
	}

	/* Check if each particle is in the box it should be in, given its
	 * coordinates and count the number of particles and check them with
	 * the total */
	for (i = 0; i < nb * nb * nb; i++)
	{
		Box *b = &world.grid[i];
		const Particle *p, *first;

		if (b->p == NULL) {
			if (b->n != 0) {
				fprintf(stderr, "Box %d: found zero particles, "
						"expected %d\n", i, b->n);
				OK = false;
			}
			continue;
		}

		first = b->p;
		p = first;
		j = 0;
		do {
			Box *correctBox = boxFromParticle(p);
			if (correctBox != b) {
				int c = (correctBox - world.grid)/sizeof(*correctBox);
				fprintf(stderr, "Particle is in box %d, "
					"should be in %ld\n", i, (correctBox -
					world.grid)/sizeof(*correctBox));
				fprintf(stderr, "numBox per dim: %d\n", nb);
				fprintf(stderr, "Pos:\t");
				fprintVector(stderr, &p->pos);
				fprintf(stderr, "\n");
				fprintf(stderr, "Actual box coords:  %d %d %d\n",
						i/nb/nb, (i/nb)%nb, i%nb);
				fprintf(stderr, "Correct box coords: %d %d %d\n",
						c/nb/nb, (c/nb)%nb, c%nb);
				OK = false;
			}
			j++;
			nParts1++;
			p = p->next;
		} while (p != first);

		if (j != b->n) {
			fprintf(stderr, "Box %d: found %d particles, "
					"expected %d\n", i, j, b->n);
			OK = false;
		}
		nParts2 += b->n;
	}

	if (nParts1 != config.numParticles)
	{
		fprintf(stderr, "1: Found a total of %d particles, "
			"should be %d\n", nParts1, config.numParticles);
		OK = false;
	}

	if (nParts2 != config.numParticles)
	{
		fprintf(stderr, "2: Found a total of %d particles, "
			"should be %d\n", nParts2, config.numParticles);
		OK = false;
	}

	return OK;
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
