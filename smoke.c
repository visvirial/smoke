/**
 * @file smoke.c
 * @author Masahiko Hyuga <mail@mhyuga.jp>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <glib.h>
#include <emmintrin.h>
#include <SDL/SDL.h>

// HEIGHT corresponds to y = 1.
// 480p
#define WIDTH (640)
#define HEIGHT (480)
// 720p
//#define WIDTH (1280)
//#define HEIGHT (720)
// 1080p
//#define WIDTH  (1920)
//#define HEIGHT (1080)

// Frames per second.
#define FPS (30)

// Background color
//#define BG_R (0x13)
//#define BG_G (0x29)
//#define BG_B (0x4D)
#define BG_R (0x00)
#define BG_G (0x00)
#define BG_B (0x00)


// How many particles create?
#define N_PARTICLES_CREATE (900)
#define MAX_PARTICLES (1000000)

static SDL_Surface *surface;

/**
 * E.O.M. is as follows,
 *   x'' = - GAMMA * x' + (N(0, ZETA), N(FY_AVG, ZETA))
 * where N(mu, sigma) is Gaussian distribution.
 */
#define GAMMA  (0.01)
#define ZETA   (0.1)
#define FY_AVG (0.005)
// Initial vy.
#define VY_INIT (0.1)
// Decay time of a particle.
#define DECAY_TIME (10.)
// Particles are created at (DX*sin(OMEGA*t) + N(0, SIGMA_X), 0).
#define DX      (0.01)
#define OMEGA   (1)
#define SIGMA_X (0.03)
// Density of a particle.
#define DENSITY (0.01)
// Typical size of a particle.
#define SIZE (0.05)

static bool draw_simple = false;

static int dens_map_s = -1;
static double **dens_map;

typedef struct{
	float r;
	float g;
	float b;
} color;

typedef struct{
	// Position vector.
	__v2df r;
	// Velocity vector.
	__v2df v;
	// Color.
	__v4sf col;
	// Lifetime^-1.
	double ltinv;
	// Left time of life.
	double t;
} particle;

static GArray *parr;

color hsv2rgb(float h, float s, float v){
	for(;h>360; h-=360);
	int hi = (int)floor(h / 60) % 6;
	float f = h / 60 - hi;
	float p = v * (1 - s);
	float q = v * (1 - f*s);
	float t = v * (1 - (1-f)*s);
	color c;
	switch(hi){
		case 0:
			c.r = v; c.g = t; c.b = p;
			break;
		case 1:
			c.r = q; c.g = v; c.b = p;
			break;
		case 2:
			c.r = p; c.g = v; c.b = t;
			break;
		case 3:
			c.r = p; c.g = q; c.b = v;
			break;
		case 4:
			c.r = t; c.g = p; c.b = v;
			break;
		case 5:
			c.r = v; c.g = p; c.b = q;
			break;
	}
	return c;
}

static void evolute(double t){
	// Create new particles.
	if(parr->len > MAX_PARTICLES){
		fprintf(stderr, "E: too many particles has been created!\n");
		exit(1);
	}
	for(int i=0; i<N_PARTICLES_CREATE/FPS; i++){
		__v2df r = {DX*sin(OMEGA*t)+SIGMA_X*sqrt(-2*log(g_random_double()))*cos(2*M_PI*g_random_double()), -SIZE};
		__v2df v = {0, VY_INIT};
		double x = g_random_double();
		double lt = -4. * DECAY_TIME * x * log(x);
		color col = hsv2rgb(36*t, 0.3, 1.);
		particle p = {r, v, {0.0, col.r, col.g, col.b}, 1./lt, lt};
		g_array_append_val(parr, p);
	}
	// Evolute with Symplectic Euler
	double dt = 1./FPS;
	__v2df vdt = {dt, dt};
	for(int i=0; i<(int)parr->len; i++){
		particle *p = &g_array_index(parr, particle, i);
		double a = g_random_double();
		double b = g_random_double();
		double t = ZETA * sqrt(-2.*log(a));
		__v2df gam = {-GAMMA, -GAMMA};
		gam *= p->v;
		__v2df fluct = {t*cos(2.*M_PI*b), t*sin(2*M_PI*b)+FY_AVG};
		__v2df f = gam + fluct;
		p->v += f * vdt;
		p->r += p->v * vdt;
		p->t -= dt;
	}
	// Delete out-of-range particles.
	for(int i=(int)parr->len-1; i>=0; i--){
		particle *p = &g_array_index(parr, particle, i);
		if(p->t < 0 || fabs(p->r[0]) > 1.1*WIDTH/HEIGHT || p->r[1] < -1.1 || p->r[1] > 2.1)
			g_array_remove_index(parr, i);
	}
	printf("\rN: %d", parr->len);
	fflush(stdout);
}

static __v4sf trans[WIDTH][HEIGHT];
static void draw(){
	uint8_t *px = (uint8_t*)surface->pixels;
	for(int x=0; x<WIDTH; x++)
		for(int y=0; y<HEIGHT; y++){
			__v4sf t = {1, 1, 1, 1};
			trans[x][y] = t;
		}
	// Draw particles
	__v4sf one = {1, 1, 1, 1};
	for(int i=0; i<(int)parr->len; i++){
		particle *p = &g_array_index(parr, particle, i);
		int x = round(.5*WIDTH + p->r[0]*HEIGHT);
		int y = round(HEIGHT*(1-p->r[1]));
		float dec = p->ltinv * p->t;
		for(int xx=-dens_map_s; xx<=dens_map_s; xx++)
			for(int yy=-dens_map_s; yy<=dens_map_s; yy++){
				if(x+xx < 0 || x+xx >= WIDTH || y+yy < 0 || y+yy >= HEIGHT) continue;
				float t = dec * dens_map[xx+dens_map_s][yy+dens_map_s];
				__v4sf vt = {t, t, t, t};
				trans[x+xx][y+yy] *= one - p->col * vt;
			}
	}
	for(int x=0; x<WIDTH; x++)
		for(int y=0; y<HEIGHT; y++){
			px[4*(x+y*WIDTH)+0] = BG_B * trans[x][y][3] + 0xFF * (1-trans[x][y][3]);
			px[4*(x+y*WIDTH)+1] = BG_G * trans[x][y][2] + 0xFF * (1-trans[x][y][2]);
			px[4*(x+y*WIDTH)+2] = BG_R * trans[x][y][1] + 0xFF * (1-trans[x][y][1]);
		}
}

static void draw_s(){
	uint8_t *px = (uint8_t*)surface->pixels;
	// Clear background.
	for(int x=0; x<WIDTH; x++)
		for(int y=0; y<HEIGHT; y++){
			px[4*(x+y*WIDTH)+0] = BG_B;
			px[4*(x+y*WIDTH)+1] = BG_G;
			px[4*(x+y*WIDTH)+2] = BG_R;
		}
	// Draw particles.
	for(int i=0; i<(int)parr->len; i++){
		particle *p = &g_array_index(parr, particle, i);
		int x = (int)(WIDTH / 2. + p->r[0] * HEIGHT);
		int y = (int)((1 - p->r[1]) * HEIGHT);
		if(x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) continue;
		float dec = p->ltinv * p->t;
		px[4*(x+y*WIDTH)+0] = BG_B * (1-dec) + 0xFF * p->col[3] * dec;
		px[4*(x+y*WIDTH)+1] = BG_G * (1-dec) + 0xFF * p->col[2] * dec;
		px[4*(x+y*WIDTH)+2] = BG_R * (1-dec) + 0xFF * p->col[1] * dec;
	}
}

static void main_loop(){
	Uint32 begin = SDL_GetTicks();
	Uint32 fps_ticks = SDL_GetTicks();
	for(long frames=0; ; frames++){
		// Process events.
		SDL_Event event;
		for(;SDL_PollEvent(&event);){
			switch(event.type){
				case SDL_KEYUP:
					if(event.key.keysym.sym == SDLK_s)
						draw_simple = !draw_simple;
					break;
				case SDL_QUIT:
					return;
			}
		}
		// Evolute.
		evolute(1.*frames/FPS);
		// Draw Window.
		if(draw_simple)
			draw_s();
		else
			draw();
		SDL_Flip(surface);
		// Calculate FPS.
		if(frames%FPS == 0){
			double fps = FPS * 1e3 / (SDL_GetTicks() - fps_ticks);
			char cap[64];
			sprintf(cap, "Smoke (FPS:%.2f)", fps);
			SDL_WM_SetCaption(cap, NULL);
			fps_ticks = SDL_GetTicks();
		}
		// Sleep if necessarry.
		Uint32 elapsed = SDL_GetTicks() - begin;
		if(0 < elapsed && elapsed*FPS < 1000)
			SDL_Delay(1000/FPS - elapsed);
		begin = SDL_GetTicks();
	}
}

//int main(int argc, char *argv[]){
int main(){
	// Initialize SDL.
	if(SDL_Init(SDL_INIT_VIDEO)){
		fprintf(stderr, "E: failed to initialize SDL library: %s\n", SDL_GetError());
		exit(1);
	}
	// Set window title.
	SDL_WM_SetCaption("Smoke (FPS:0.00)", NULL);
	// Initialize surface.
	surface = SDL_SetVideoMode(WIDTH, HEIGHT, 32, SDL_HWSURFACE | SDL_DOUBLEBUF);
	if(surface == NULL){
		fprintf(stderr, "E: failed to initialize SDL surface: %s\n", SDL_GetError());
		exit(1);
	}
	
	// Initialize.
	parr = g_array_new(FALSE, FALSE, sizeof(particle));
	dens_map_s = ceil(SIZE * HEIGHT);
	dens_map = malloc((2*dens_map_s+1)*sizeof(double*));
	for(int i=0; i<2*dens_map_s+1; i++)
		dens_map[i] = malloc((2*dens_map_s+1)*sizeof(double));
	for(int x=-dens_map_s; x<=dens_map_s; x++){
		for(int y=-dens_map_s; y<=dens_map_s; y++){
			double q = 3.;
			double r2 = (double)(x*x+y*y) / dens_map_s / dens_map_s;
			double a = DENSITY / (1. - exp(-q));
			double b = DENSITY - a;
			dens_map[x+dens_map_s][y+dens_map_s] = (r2>1 ? 0 : a*exp(-q*r2)+b);
		}
	}
	
	// Call main loop.
	main_loop();
	printf("\n");
	
	// Cleanup.
	g_array_free(parr, TRUE);
	for(int i=0; i<2*dens_map_s+1; i++)
		free(dens_map[i]);
	free(dens_map);
	SDL_Quit();
	return 0;
}

