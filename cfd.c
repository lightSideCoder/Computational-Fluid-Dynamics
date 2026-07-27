#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>
#include <string.h>
#include <time.h>

#define WIDTH 800
#define HEIGHT 800
#define DT (1.0f / 60.0f)
#define DIV 10
#define COLS ((WIDTH) / (DIV) + 2)
#define ROWS ((HEIGHT) / (DIV) + 2)
#define CELLS ((COLS) * (ROWS))
#define DENS 1
#define ITER 60
#define OVERR 1.9
#define U_FIELD 1
#define V_FIELD 2
#define S_FIELD 3
#define G 0//9.81
#define NUM_P 20000


typedef struct {
	float x, y;
} Particle;


float U[CELLS] = {0};
float V[CELLS] = {0};
float newU[CELLS] = {0};
float newV[CELLS] = {0};
float P[CELLS] = {0};
float S[CELLS] = {0};
float M[CELLS] = {0};
float newM[CELLS] = {0};

Particle particles[NUM_P];


void integrate();
void solveIncompressibility();
void extrapolate();
float sampleField(float x, float y, int field, float *M);
float avgU(int i, int j);
float avgV(int i, int j);
void advectVelocity(float *M);
void advectSmoke(float *M);
void simulate(float *M);
void mouse();


int main() {
	srand(time(NULL));

	for (int i = 0; i < NUM_P; i++) {
		particles[i].x = rand() % WIDTH;
		particles[i].y = rand() % HEIGHT;
	}

//----------------------------------------------------------------------------------------------------------------------
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window *pwindow = SDL_CreateWindow("CFD", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
	SDL_Renderer *renderer = SDL_CreateRenderer(pwindow, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
	int running = 1;
//----------------------------------------------------------------------------------------------------------------------

	while(running) {
		SDL_Event event;
		while(SDL_PollEvent(&event)) {
			switch(event.type) {
				case SDL_QUIT:
					running = 0;
					break;
			}
		}
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
		SDL_RenderClear(renderer);

		SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);

		mouse();
		simulate(M);

		for (int i = 0; i < NUM_P; i++) {
			float u = sampleField(particles[i].x, particles[i].y, U_FIELD, M);
			float v = sampleField(particles[i].x, particles[i].y, V_FIELD, M);
			particles[i].x += u * DT;
			particles[i].y += v * DT;

			SDL_RenderDrawPointF(renderer, particles[i].x, particles[i].y);
		}
		
		SDL_RenderPresent(renderer);
	}

	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(pwindow);
	SDL_Quit();

	return 0;
}

void integrate() {
	int n = ROWS;
	for (int i = 1; i < COLS; i++) {
	for (int j = 1; j < n-1; j++) {
		if (S[i*n+j] != 0 && S[i*n +j-1] != 0) {
			V[i*n+j] += G * DT;
		}
	}
	}
}

void solveIncompressibility() {
	int n = ROWS;
	float cp = DENS * DIV / DT;
	for (int iter = 0; iter < ITER; iter++) {
		for (int i = 1; i < COLS-1; i++) {
		for (int j = 1; j < ROWS-1; j++) {
			if (S[i*n+j] == 0) continue;
			
			//If the following points are walls or objects their value is 0 if they flow freely it's 1
			float sx0 = S[(i-1)*n+j];
			float sx1 = S[(i+1)*n+j];
			float sy0 = S[i*n+j-1];
			float sy1 = S[i*n+j+1];
			float s = sx0 + sx1 + sy0 + sy1; //Get the value we need to divide the divergence by, that we then add to / substract from the velocity vectors
			if (s == 0)  continue;

			//Divergence is how much flows in to / out of the cell. Incompressibility means this value is 0
			float divergence = U[(i+1)*n+j] - U[i*n+j] + V[i*n+j+1] - V[i*n+j];

			float p = -divergence / s; //The amount that flows in to / out of the cell gets divided by the number of vectors surrounding it, that aren't walls
			p *= OVERR;
			P[i*n+j] += cp * p;

			U[i*n+j] -= sx0 * p;	//If the point is actually flowing and not a wall we add / subtract p to/from it
			U[(i+1)*n+j] += sx1 * p;//Doing that with all velocities surrounding the cell, the value of divergence will get close to zero
			V[i*n+j] -= sy0 * p;	//Through iteration the divergence over the whole domain gets close to zero
			V[i*n+j+1] += sy1 * p;
		}
		}
	}
}


void extrapolate() {
	int n = ROWS;
	for (int i = 0; i < COLS; i++) {
		U[i*n] = U[i*n+1];			//We take the velocity U at the top border and make it equal to the second row
		U[i*n+n-1] = U[i*n+n-2];		//We take the velocity U at the bottom border and make it equal to the row one above it
	}
	for (int j = 0; j < ROWS; j++) {
		V[j] = V[n+j];				//We take the velocity V at the left border and make it equal to the velocity one column to the right
		V[(COLS-1)*n+j] = V[(COLS-2)*n+j];	//We take the velocity V at the right border and make it equal to the velocity one column to the left
	}
}


float sampleField(float x, float y, int field, float *M) {
	int n = ROWS;
	float h = DIV;
	float h1 = 1.0f / h;
	float h2 = h * .5;

	x = fminf(fmaxf(x, h), (COLS-1.001f) * h); //If the argument x is out of bounds (bigger than COLS*h, smaller than h), we will set x to either COLS*h or h respectively
											   //If it is in bounds just leave it as it is
	y = fminf(fmaxf(y, h), (ROWS-1.001f) * h); //If the argument y is out of bounds (bigger than ROWS*h, smaller than h), we will set < to either ROWS*h or h respectively

	float dx = 0;
	float dy = 0;

	float *F; //Initialize Pointer to prepare for our array-swap

	//Depending on the argument "field" we swap either the values from U, V, or M with F. And either set dy or dx or both to half the height of a cell, depending on the argument
	switch (field) {
		case U_FIELD: F = U; dy = h2; break;
		case V_FIELD: F = V; dx = h2; break;
		case S_FIELD: F = M; dx = h2; dy = h2; break;
	}


/*	Here starts the BI-LINEAR INTERPOLATION:

	We wanna find out:	In which cell are we? -> x0, y0
				How far are we from the cell borders --> tx, ty
				Whats the neighbouring cell borders --> x1, y1
*/
	//x0 ist basically die linke Grenze der Zelle, in der wir uns befinden
	float x0 = (floor((x-dx)*h1) < COLS-1) ? floor((x-dx)*h1) : COLS-1;
	
	//tx is the distance in px from that border
	float tx = ((x-dx) - x0*h) * h1; //x-dx ist wieder x-0 bzw x-5 (je nachdem welches field wir "samplen") davon ziehen wir x0*10 ab, das multiplizieren wir dann mit 0.1

	//x1 is the right cell border of our cell
	float x1 = (x0+1 < COLS-1) ? x0+1 : COLS-1;


	//y0 is the top border of our cell
	float y0 = (floor((y-dy)*h1) < ROWS-1) ? floor((y-dy)*h1) : ROWS-1;

	//ty is the distance from that border
	float ty = ((y-dy) - y0*h) * h1;

	//y1 is the bottom border of our cell
	float y1 = (y0+1 < ROWS-1) ? y0+1 : ROWS-1;

	float sx = 1 - tx;
	float sy = 1 - ty;

	//The closer to the corner of the cell we are, the stronger the value stored in F[] at the cell corresponding to that corner will be present in val
	float val = sx*sy * F[(int)(x0*n + y0)] +
		    tx*sy * F[(int)(x1*n + y0)] +
		    tx*ty * F[(int)(x1*n + y1)] +
		    sx*ty * F[(int)(x0*n + y1)];

	return val;
}

float avgU(int i, int j) {
	float u = (U[i*ROWS+j-1] + U[i*ROWS+j] +
		  U[(i+1)*ROWS+j-1] + U[(i+1)*ROWS+j]) * .25;
	return u;
}

float avgV(int i, int j) {
	float v = (V[(i-1)*ROWS+j] + V[i*ROWS+j] +
		  V[(i-1)*ROWS+j+1] + V[i*ROWS+j+1]) * .25;
	return v;
}

void advectVelocity(float *M) {
	memcpy(newU, U, CELLS*sizeof(float));
	memcpy(newV, V, CELLS*sizeof(float));

	int n = ROWS;
	int h = DIV;
	int h2 = h * .5;

	for (int i = 1; i < COLS; i++) {
	for (int j = 1; j < ROWS; j++) {
		if (S[i*n+j] != 0 && S[(i-1)*n+j] != 0 && j < n -1) {
			float x = i * h;	//10, 20, 30, ..., 800
			float y = j * h + h2;	//15, 25, ..., 805
			float u = U[i*n+j];
			float v = avgV(i, j);
			x -= u * DT;
			y -= v * DT;
			u = sampleField(x, y, U_FIELD, M);
			newU[i*n+j] = u;
		}
		if (S[i*n+j] != 0 && S[i*n+j-1] != 0 && i < COLS -1) {
			float x = i * h + h2;
			float y = j * h;
			float u = avgU(i, j);
			float v = V[i*n+j];
			x -= u * DT;
			y -= v * DT;
			v = sampleField(x, y, V_FIELD, M);
			newV[i*n+j] = v;

		}
	}
	}
	memcpy(U, newU, CELLS*sizeof(float)); //Array kopieren
	memcpy(V, newV, CELLS*sizeof(float)); //Array kopieren
}

void advectSmoke(float *M) {

	memcpy(newM, M, CELLS*sizeof(float)); //Array kopieren
	
	int n = ROWS;
	int h = DIV;
	int h2 = DIV * .5;

	for (int i = 1; i < COLS-1; i++) {
	for (int j = 1; j < n-1; j++) {
		if (S[i*n+j] != 0) {
			float u = (U[i*n+j] + U[(i+1)*n+j]) * .5;
			float v = (V[i*n+j] + V[i*n+j+1]) * .5;
			float x = i*h + h2 - u * DT;
			float y = j*h + h2 - v * DT;

			newM[i*n+j] = sampleField(x, y, S_FIELD, M);
		}
	}
	}
	memcpy(M, newM, CELLS*sizeof(float)); //Array kopieren
}

void simulate(float *M) {
	for (int i = 0; i < COLS; i++) {
	for (int j = 0; j < ROWS; j++) {
		if (i == 0 || j == 0 || i == COLS-1 || j == ROWS-1) {
			S[i*ROWS+j] = 0;
		} else S[i*ROWS+j] = 1;
	}
	}

	integrate();

	for (int i = 0; i < CELLS; i++) {
		P[i] = 0;
	}

	solveIncompressibility();
	extrapolate();
	advectVelocity(M);
	advectSmoke(M);
}

void mouse(void)
{
    static int lastX = -1;
    static int lastY = -1;

    int mx, my;
    Uint32 buttons = SDL_GetMouseState(&mx, &my);

    if (!(buttons & SDL_BUTTON(SDL_BUTTON_LEFT))) {
        lastX = mx;
        lastY = my;
        return;
    }

    // mouse velocity in pixels/frame
    float dx = mx - lastX;
    float dy = my - lastY;

    lastX = mx;
    lastY = my;

    // convert to grid coordinates
    int i = mx / DIV;
    int j = my / DIV;

    for (int di = -2; di <= 2; di++) {
        for (int dj = -2; dj <= 2; dj++) {

            int ii = i + di;
            int jj = j + dj;

            if (ii < 1 || ii >= COLS-2 || jj < 1 || jj >= ROWS-2) {
                continue;
			}

            // inject smoke
	    	M[ii*ROWS + jj] = 1.0f;

            // inject momentum
            U[ii*ROWS + jj] += dx * 5.0f;
            U[(ii+1)*ROWS + jj] += dx * 5.0f;

            V[ii*ROWS + jj] += dy * 5.0f;
            V[ii*ROWS + jj + 1] += dy * 5.0f;
        }
    }
}
