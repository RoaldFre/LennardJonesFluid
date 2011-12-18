#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "vmath.h"
#include "system.h"
#include "measure.h"

#define NOT_USED(x) ( (void)(x) )

static void collideWalls(int ix, int iy, int iz);
static void reboxParticles(void);
static Box  *boxFromParticle(const Particle *p);
static Box  *boxFromIndex(int nx, int ny, int nz);
static Box  *boxFromNonPeriodicIndex(int nx, int ny, int nz);
static void addToBox(Particle *p, Box *b);
static void removeFromBox(Particle *p, Box *from);

static void verlet(void);
static void thermostat(void);
static void calculateAcceleration(void);
static void lennardJonesForce(Particle *p1, Particle *p2);


static double randNorm(void);




/* GLOBALS */

World world;
Config config;

double sim_time = 0;


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
 * 'config.thermostatTemp'. Velocities are then shifted to set the total 
 * momentum to zero. */
void fillWorld()
{
	Particle *ps = world.parts;
	double worldSize = config.numBox * config.boxSize;
	double stdev = sqrt(config.initialTemp);
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

/* Execute a given function for every discinct pair of particles that are 
 * within the same box, or in adjacent boxes..
 * Arguments:
 *  - Function pointer to function that will be fed all the particle pairs.
 *  - Pointer to data that will be supplied to said function.
 */
void forEveryPair(void (*f)(Particle *p1, Particle *p2, void *data), void *data)
{
	int nb = config.numBox;
	int n1, n2;

	assert(sanityCheck());

	if (nb < 3) {
		/* Brute force. Reason: see comment below */
		for (int i = 0; i < config.numParticles; i++) {
			Particle *p1 = &world.parts[i];
			for (int j = i + 1; j < config.numParticles; j++) {
				Particle *p2 = &world.parts[j];
				(*f)(p1, p2, data);
			}
		}
		return;
	}


	/* Loop over all boxes */
	for (int ix = 0; ix < nb; ix++)
	for (int iy = 0; iy < nb; iy++)
	for (int iz = 0; iz < nb; iz++) {
		/* Loop over all particles in this box */
		Box *box = boxFromIndex(ix, iy, iz);
		Particle *p = box->p;

		/* Loop over every partner of the i'th particle 'p' from the 
		 * box 'box' */
		n1 = box->n;
		for (int i = 0; i < n1; i++) { /* i'th particle in box */
			Particle *p2 = p->next;
			for (int j = i + 1; j < n1; j++) {
				(*f)(p, p2, data);
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
				if (b <= box)
					continue;
					/* if b == box: it's our own box!
					 * else: only check boxes that have 
					 * a strictly larger pointer value 
					 * to avoid double work. */
				p2 = b->p;
				n2 = b->n;
				for (int j = 0; j < n2; j++) {
					(*f)(p, p2, data);
					p2 = p2->next;
				}
				assert(p2 == b->p);
			}

			p = p->next;
		}

		assert(p == box->p);
	}
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
	for (int i = 0; i < config.numParticles; i++) {
		Particle *p = &world.parts[i];
		Vec3 tmp;

		/* vel(t + dt) = vel(t + dt/2) + acc(t + dt)*dt/2 */
		scale(&p->acc, dt/2, &tmp);
		add(&p->vel, &tmp, &p->vel);
	}
}

static void accelerationWorker(Particle *p1, Particle *p2, void *data)
{
	NOT_USED(data);
	lennardJonesForce(p1, p2);
}

static void calculateAcceleration()
{
	/* Reset acceleration */
	for (int i = 0; i < config.numParticles; i++) {
		world.parts[i].acc.x = 0;
		world.parts[i].acc.y = 0;
		world.parts[i].acc.z = 0;
	}

	forEveryPair(&accelerationWorker, NULL);
}

/* Returns the shortest vector that points from v1 to v2, taking into 
 * account the periodic boundary conditions. The particle MUST be within 
 * the correct bounds. */
Vec3 nearestImageVector(Vec3 *v1, Vec3 *v2)
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

	//assert(fabs(-24*(2/dr12 - 1/dr6)) < 1e10);

	scale(&drVec, -24*(2/dr14 - 1/dr8), &Fi);
	add(&p1->acc, &Fi, &p1->acc);
	sub(&p2->acc, &Fi, &p2->acc);
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
	double T0  = config.thermostatTemp;
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
	sim_time += config.timeStep;
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

