#include <stdio.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <math.h>
#include "main.h"
#include "render.h"
#include "system.h"

#define SPHERE_SLICES 8

#define SCREEN_W 500
#define SCREEN_H 500

typedef struct {
	GLfloat x, y, z;
} Vertex3;

const GLfloat light_pos[]  = {2.0, 1.0, 2.0, 0.0};
const GLfloat light_diff[] = {1.0, 1.0, 1.0, 0.0};
const GLfloat light_spec[] = {1.0, 0.0, 0.0, 0.0};
const GLfloat light_ambi[] = {0.8, 0.8, 0.8, 0.0};

static int numVertices;
static int numIndices;
static Vertex3 *sphereVertex;
static GLushort *sphereIndex;

static SDL_Surface *surface;
static double angle;

static void createSphere(int slices, int *numVert, Vertex3 **vertices, int *numInd,
		GLushort **indices);
static void gluPerspective(GLfloat fovy, GLfloat aspect, GLfloat zNear, 
		GLfloat zFar);
static void calcFps(void);

static void calcFps()
{
	static int tock = 0, frames;
	int tick;
	char string[32];

	frames++;
	tick = SDL_GetTicks();

	if (tick - tock > 1000)
	{
		tock = tick;
		sprintf(string, "%u FPS", frames);
		SDL_WM_SetCaption(string, string);
		frames = 0;
	}
	return;
}

/* Update the rendered image and handle events
 * Return: true if everything went fine, false if user requested to quit. */
bool stepGraphics()
{
	SDL_Event event;
	while (SDL_PollEvent(&event))
	{
		if (event.type == SDL_QUIT) {
			printf("\nRequested Quit.\n\n");
			return false;
		}
		else if (event.type == SDL_KEYDOWN)
		{
			switch (event.key.keysym.sym)
			{
			case SDLK_ESCAPE:
				printf("\nRequested Quit.\n\n");
				return false;
				break;
			case SDLK_LEFT:
				angle--;
				break;
			case SDLK_RIGHT:
				angle++;
				break;
			case SDLK_UP:
				config.timeStep *= 1.2;
				printf("Time step: %f\n", config.timeStep);
				break;
			case SDLK_DOWN:
				config.timeStep /= 1.2;
				printf("Time step: %f\n", config.timeStep);
				break;
			case SDLK_SPACE:
				config.thermostatTemp *= 1.2;
				printf("Temperature: %f\n", config.thermostatTemp);
				break;
			case SDLK_BACKSPACE:
				config.thermostatTemp /= 1.2;
				printf("Temperature: %f\n", config.thermostatTemp);
				break;
			case SDLK_RETURN:
				SDL_WM_ToggleFullScreen(surface);
				break;
			default:
				break;
			}
		}
	}

	render();
	calcFps();

	return true;
}


static void gluPerspective(GLfloat fovy, GLfloat aspect, GLfloat zNear, 
		GLfloat zFar)
{
	GLfloat xMin, xMax, yMin, yMax;

	yMax = zNear * tan(fovy * M_PI / 360.0);
	yMin = -yMax;

	xMin = yMin * aspect;
	xMax = yMax * aspect;

	glFrustum(xMin, xMax, yMin, yMax, zNear, zFar);
}

int initRender(void)
{
	int flags = 0;
	double s = config.numBox * config.boxSize;
	const SDL_VideoInfo *vidinfo;

	if (SDL_Init(SDL_INIT_VIDEO) < 0)
		die(SDL_GetError());
	
	vidinfo = SDL_GetVideoInfo();
	if (vidinfo == NULL)
		die(SDL_GetError());
	
	flags |= SDL_OPENGL;
	flags |= SDL_HWPALETTE;
	flags |= (vidinfo->hw_available ? SDL_HWSURFACE : SDL_SWSURFACE);
	if (vidinfo->blit_hw) flags |= SDL_HWACCEL;

	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 0);
	SDL_GL_SetAttribute(SDL_GL_SWAP_CONTROL, 1);

	surface = SDL_SetVideoMode(SCREEN_W, SCREEN_H, 24, flags);
	if (surface == NULL)
		die(SDL_GetError());

	SDL_EnableKeyRepeat(1, 
			SDL_DEFAULT_REPEAT_INTERVAL);

	/*SDL_WM_ToggleFullScreen(surface);	*/

	atexit(SDL_Quit);

	/* OpenGL Init */
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);

	glEnable(GL_LIGHTING);

	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diff);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_spec);
	glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambi);

	glClearColor(1.0, 1.0, 1.0, 0.0);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

	createSphere(SPHERE_SLICES, &numVertices, &sphereVertex, &numIndices, &sphereIndex);
	glVertexPointer(3, GL_FLOAT, sizeof(Vertex3), sphereVertex);
	glNormalPointer(   GL_FLOAT, sizeof(Vertex3), sphereVertex);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(35, SCREEN_W/(double)SCREEN_H, s/2, 4*s);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	return 0;
}

int render(void)
{
	int i;
	double s;

	s = config.numBox * config.boxSize;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

/*	angle += 0.01;*/

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0, 0, -s*2.5);
	glRotatef(angle, 0, 1, 0);

	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINE_LOOP);
		glVertex3f(-s/2, -s/2, -s/2);
		glVertex3f(-s/2, -s/2, +s/2);
		glVertex3f(-s/2, +s/2, +s/2);
		glVertex3f(-s/2, +s/2, -s/2);
	glEnd();

	glBegin(GL_LINE_LOOP);
		glVertex3f(+s/2, -s/2, -s/2);
		glVertex3f(+s/2, -s/2, +s/2);
		glVertex3f(+s/2, +s/2, +s/2);
		glVertex3f(+s/2, +s/2, -s/2);
	glEnd();

	glColor3f(0.0, 0.7, 0.0);
	glPointSize(2);
	for (i = 0; i < config.numParticles; i++)
	{
		Particle *p = &world.parts[i];

		glPushMatrix();
			glTranslatef(p->pos.x - s/2, p->pos.y - s/2, p->pos.z - s/2);
			glScalef(config.radius, config.radius, config.radius);
			glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_SHORT, sphereIndex);
		glPopMatrix();
	}

	SDL_GL_SwapBuffers();

	return 0;
}

static void createSphere(int slices, int *numVert, Vertex3 **vertices, int *numInd,
		GLushort **indices)
{
	int i, j, k;
	double x, y, z;
	double r;
	int stacks;
	Vertex3 *vert;
	GLushort *ind;

	stacks = slices;
	slices *= 2;

	/* Plus two for the poles */
	*numVert = (stacks - 1) * slices + 2;
	*vertices = calloc(*numVert, sizeof(Vertex3));
	vert = *vertices;

	/* All but the top and bottom stack */
	for (i = 1; i < stacks; i++)
	{
		double phi = M_PI * i / (double) stacks - 2*M_PI;
		
		z = cos(phi);
		r = sqrt(1 - z*z);

		for (j = 0; j < slices; j++)
		{
			double theta = 2*M_PI*j/(double) slices;
			x = r * sin(theta);
			y = r * cos(theta);

			vert[(i-1) * slices + j + 1].x = x;
			vert[(i-1) * slices + j + 1].y = y;
			vert[(i-1) * slices + j + 1].z = z;
		}
	}

	/* Top and bottom */
	vert[0].x = 0;
	vert[0].y = 0;
	vert[0].z = 1;

	vert[*numVert-1].x = 0;
	vert[*numVert-1].y = 0;
	vert[*numVert-1].z = -1;

	*numInd = (stacks - 1) * slices * 6;
	*indices = calloc(*numInd, sizeof(GLushort));
	ind = *indices;

	k = 0;

	for (i = 1; i < slices; i++)
	{
		ind[k++] = 0;
		ind[k++] = i;
		ind[k++] = i+1;
	}
	ind[k++] = 0;
	ind[k++] = 1;
	ind[k++] = slices;
	
	for (i = 0; i < slices - 1; i++)
	{
		ind[k++] = *numVert - 1;
		ind[k++] = (*numVert - 1 - slices) + i;
		ind[k++] = (*numVert - 1 - slices) + i + 1;
	}
	ind[k++] = *numVert - 1;
	ind[k++] = *numVert - 1 - 1;
	ind[k++] = *numVert - 1 - slices + 0;

	for (i = 1; i < stacks - 1; i++)
	{
		int base = 1 + (i - 1) * slices;

		for (j = 0; j < slices - 1; j++)
		{
			ind[k++] = base + j;
			ind[k++] = base + slices + j;
			ind[k++] = base + slices + j + 1;

			ind[k++] = base + j;
			ind[k++] = base + j + 1;
			ind[k++] = base + slices + j + 1;
		}

		ind[k++] = base;
		ind[k++] = base + slices - 1;
		ind[k++] = base + slices;

		ind[k++] = base + slices - 1;
		ind[k++] = base + slices;
		ind[k++] = base + slices + slices - 1;
	}

	return;
}
