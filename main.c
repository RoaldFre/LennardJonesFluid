#define _GNU_SOURCE

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "main.h"
#include "vmath.h"
#include "system.h"
#include "render.h"

#define FILENAME "/tmp/data.txt"

static void parseArguments(int argc, char **argv);
static void plotHeader(FILE *stream);
static void plotFooter(FILE *stream);
static int plot(void);

/* Defaults */
#define DEF_TIMESTEP 			0.0001
#define DEF_TEMPERATURE 		2.0
#define DEF_COUPLING_TIMESTEP_FACTOR 	1000
#define DEF_LJ_TRUNCATION 		2.5
#define DEF_MEASUREMENT_INTERVAL	0.2
#define DEF_MEASUREMENT_WAIT 		10.0
#define DEF_PARTICLES_PER_RENDER 	10000
#define DEF_RENDER_RADIUS 		1.0

static void printUsage(void)
{
	printf("Usage: main <number of particles> <particle density> [flags]\n");
	printf("Flags:\n");
	printf(" -t <flt>  length of Time steps\n");
	printf("             default: %f\n", DEF_TIMESTEP);
	printf(" -T <flt>  Temperature\n");
	printf("             default: %f (Boltzmann cte = mass = 1)\n", DEF_TEMPERATURE);
	printf(" -c <flt>  thermal bath Coupling: relaxation time (zero to disable)\n");
	printf("             default: %d * timestep\n", DEF_COUPLING_TIMESTEP_FACTOR);
	printf(" -l <flt>  Lennard-Jones truncation Length\n");
	printf("             default: %f\n", DEF_LJ_TRUNCATION);
	printf(" -b <num>  number of Boxes per dimension\n");
	printf("             default: max so that boxsize >= L-J truncation length\n");
	printf(" -s <int>  accumulate <int> measurement Samples\n");
	printf("             default: Don't sample, loop forever\n");
	printf(" -i <flt>  time Interval of <flt> between measurements\n");
	printf("             default: %f\n", DEF_MEASUREMENT_INTERVAL);
	printf(" -w <flt>  Wait for a time <flt> before starting the measurements\n");
	printf("             default: %f\n", DEF_MEASUREMENT_WAIT);
	printf(" -j <int>  perform <int> physics steps between rendering frames.\n");
	printf("             default: %d/(number of particles)\n", DEF_PARTICLES_PER_RENDER);
	printf(" -r        Render\n");
	printf(" -R <flt>  Radius of the particles when rendering\n");
	printf("             default: %f\n", DEF_RENDER_RADIUS);
	printf(" -v <int>  Verbose: dump statistics every <int> iterations\n");
}

static void parseArguments(int argc, char **argv)
{
	int c;

	/* defaults */
	config.verbose = 0;
	config.measureSamples = -1; /* loop indefinitely */
	config.truncateLJ 	= DEF_LJ_TRUNCATION;
	config.timeStep 	= DEF_TIMESTEP;
	config.measureInterval  = DEF_MEASUREMENT_INTERVAL;
	config.measureWait	= DEF_MEASUREMENT_WAIT;
	config.temperature	= DEF_TEMPERATURE;
	config.radius		= DEF_RENDER_RADIUS;

	/* guards */
	config.renderSteps = -1;
	config.numBox = -1;
	config.thermostatTau = -1;

	while ((c = getopt(argc, argv, ":t:T:c:l:b:s:i:w:j:rR:v:h")) != -1)
	{
		switch (c)
		{
		case 't':
			config.timeStep = atof(optarg);
			if (config.timeStep <= 0)
				die("Invalid timestep %f\n",
						config.timeStep);
			break;
		case 'T':
			config.temperature = atof(optarg);
			if (config.temperature < 0)
				die("Invalid temperature %f\n",
						config.temperature);
			break;
		case 'c':
			config.thermostatTau = atof(optarg);
			if (config.thermostatTau <= 0)
				die("Invalid thermostat relaxation time %f\n",
						config.thermostatTau);
			break;
		case 'l':
			config.truncateLJ = atof(optarg);
			if (config.truncateLJ <= 0)
				die("Invalid LJ truncation length %f\n",
						config.truncateLJ);
			break;
		case 'b':
			config.numBox = atoi(optarg);
			if (config.numBox <= 0)
				die("Invalid number of boxes %d\n",
						config.numBox);
			break;
		case 's':
			config.measureSamples = atol(optarg);
			if (config.measureSamples < 0)
				die("Invalid number of samples %d\n",
						config.measureSamples);
			break;
		case 'i':
			config.measureInterval = atof(optarg);
			if (config.measureInterval < 0)
				die("Invalid interval time %f\n", 
						config.measureInterval);
			break;
		case 'w':
			config.measureWait = atof(optarg);
			if (config.measureWait < 0)
				die("Invalid wait time %f\n", config.measureWait);
			break;
		case 'j':
			config.renderSteps = atol(optarg);
			if (config.renderSteps < 0)
				die("Invalid number of renderer steps %d\n",
						config.renderSteps);
			break;
		case 'r':
			config.render = true;
			break;
		case 'R':
			config.radius = atof(optarg);
			if (config.radius <= 0)
				die("Invalid radius %f\n", config.radius);
			break;
		case 'v':
			config.verbose = atoi(optarg);
			if (config.verbose <= 0)
				die("Verbose: invalid number of iterations %d\n",
						config.verbose);
			break;
		case 'h':
			printUsage();
			exit(0);
			break;
		case ':':
			printUsage();
			die("Option -%c requires an argument\n", optopt);
			break;
		case '?':
			printUsage();
			die("Option -%c not recognized\n", optopt);
			break;
		default:
			/* XXX */
			break;
		}
	}

	argc -= optind;
	argv += optind;

	if (argc < 2) {
		printUsage();
		die("Not enough required arguments!\n");
	}

	config.numParticles = atoi(argv[0]);
	double density = atof(argv[1]);

	if (config.measureInterval < config.timeStep)
		die("The interval time between measurements %f is smaller "
			"than the time step %f!\n", config.measureInterval, 
						    config.timeStep);

	double worldSize = cbrt(config.numParticles / density);

	if (config.numBox == -1)
		config.numBox = worldSize / config.truncateLJ;

	config.boxSize = worldSize / config.numBox;

	if (config.boxSize < config.truncateLJ)
		die("The boxsize (%f) is smaller than the L-J truncation "
			"radius (%f)!\n", config.boxSize, config.truncateLJ);

	if (config.thermostatTau < 0)
		config.thermostatTau = 1000 * config.timeStep;

	if (config.renderSteps < 0)
		config.renderSteps = 1 + 10000 / config.numParticles;

	printf("Boxsize: %f\n", config.boxSize);
}

void die(const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	exit(1);
}

/* Advance the simulation by one time step. Render if neccesary. Return 
 * false if the user wants to quit. */
static bool stepSimulation(void) {
	static long stepsSinceRender = 0;
	static int  stepsSinceVerbose = 0;

	stepWorld();

	if (config.verbose > 0) {
		stepsSinceVerbose++;
		if (stepsSinceVerbose > config.verbose) {
			stepsSinceVerbose = 0;
			dumpStats();
		}
	}

	if (config.render) {
		stepsSinceRender++;
		if (stepsSinceRender > config.renderSteps) {
			stepsSinceRender = 0;
			return stepGraphics();
		}
	}

	return true;
}

int main(int argc, char **argv)
{
	bool keepGoing = true;

	parseArguments(argc, argv);

	allocWorld();
	fillWorld();
	
	if (config.render)
		initRender();

	if (config.measureSamples < 0) {
		/* Loop forever, or until the user quits the renderer */
		while (stepSimulation());
	} else {
		printf("Waiting for system to relax.\n");
		for (double t = 0; keepGoing && t < config.measureWait; t += config.timeStep) {
			keepGoing = stepSimulation();
			if (fmod(t, config.measureWait / 100) < config.timeStep) {
				printf("\rRelax time %13f of %f", t + config.measureWait/100, 
								  config.measureWait);
				fflush(stdout);
			}
		}

		/* Perform the measurements */
		printf("\nStarting measurement.\n");
		FILE *outfile = fopen(FILENAME, "w");
		//plotHeader(outfile);
		double intervalTime = 0;
		for (long sample = 0; keepGoing && sample < config.measureSamples; sample++) {
			while (keepGoing && intervalTime <= config.measureInterval) {
				keepGoing = stepSimulation();
				intervalTime += config.timeStep;
			}
			if (!keepGoing)
				break;

			/* Check for numerical drift (or bugs) before 
			 * commiting measurement. */
			if (!physicsCheck()) {
				fprintf(stderr, "Something went wrong! Dumping the "
						"measurement I have till now.\n");
				dumpMeasurement(outfile);
				die("You broke physics!\n");
			}

			printf("\rMeasured sample %ld/%ld", sample + 1, config.measureSamples);
			fflush(stdout);
			sampleMeasurement();
			intervalTime -= config.measureInterval;
		}
		printf("\n");
		
		dumpMeasurement(outfile);
		//plotFooter(outfile);
		fclose(outfile);
		//plot();
	}

	freeWorld();
	return 0;
}

static void plotHeader(FILE *stream)
{
	fprintf(stream, "data = ");
}
static void plotFooter(FILE *stream)
{
	fprintf(stream, "\n");
	fprintf(stream, "plot(data(:,1), data(:,2))\n");
}
static int plot()
{
#if 0
	if (fork() > 0)
		return; /* The parent does nothing  */

	/* We are the child */
	setsid(); /* Detatch from terminal */
	/* Close inherited file descriptors */
	fclose(stdin);
	fclose(stdout);
	fclose(stderr);
#endif

	return system("octave -q --persist "FILENAME);
}
