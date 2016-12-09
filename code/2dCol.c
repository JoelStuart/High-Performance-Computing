/**************************
* Two dimensional elastic collisions - Project 11
* 
* Author:  Joel Stuart - 43203714 - 19/08/2016
* __________________
* Use guide:
* Modify define values to change number of objects or timesteps 
* or width of target cloud. NOTE: Whiteball / bullet object is 
* set in code - may have to go digging.
*
* Program outputs to static filename defined below. Read this output
* for verification
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Radius of objects
static double rad = 3;
//Filename to output to.
static char* filename = "output.txt";

//These values may be modified
#define NUM_OBJ 25
#define NUM_STEPS 10
#define X_MAX 40
#define Y_MAX 30
#define X_VEC_MAX 5
#define Y_VEC_MAX 4

/**An obj contains a array of doubles of size 4
   index
   0 - x position
   1 - y position
   2 - x velocity
   3 - y velocity
*/
typedef struct {
	double* pos;
} obj;

static void time_step(obj** stepCollec);
static obj* handle_step(obj* in, int s);
obj* gen_first_step();
int rand_int(int n);
double cpu_time(void);
void print_step(obj* step, int d);
double* collide(double* i, double* j);
void print_collision(int a, int b, int d);
void print_timestep(int d);

void main(void) {

	printf("Started up...\n");
	double start = cpu_time();
	remove(filename);
	//Generate the objects in the initial time step
	obj* firstStep = gen_first_step();
	
	//Allocate memory for the collection of time steps
	obj** stepCollec = malloc(sizeof(obj[NUM_OBJ])*NUM_STEPS);
	stepCollec[0] = firstStep;

	//Perform time stepping
	time_step(stepCollec);
	
	double stop = cpu_time()-start;
	
	printf("Finished in %f seconds. Exiting...\n", stop);


  return;
}

void time_step(obj** stepCollec) {
	int i;
	for(i=1; i < NUM_STEPS+1; i++){
		stepCollec[i] = handle_step(stepCollec[i-1], i);
		//free(stepCollec[i-1]);
	}
	
	return;
	
}

obj* handle_step(obj* in, int s) {
	int i,j;
	print_timestep(s);
	//Do time advancement on objects
	obj* objList = malloc(sizeof(obj)*2*NUM_OBJ);
	for(i=0; i < NUM_OBJ; i++){
		double tAx = in[i].pos[0] + in[i].pos[2];
		double tAy = in[i].pos[1] + in[i].pos[3];
		double* temp = malloc(sizeof(double)*4);
		temp[0] = tAx;
		temp[1] = tAy;
		temp[2] = in[i].pos[2];
		temp[3] = in[i].pos[3];
		obj a;
		a.pos = temp;
		objList[i] = a;

	}

	//Check collisions between every object
	for(i=0; i < NUM_OBJ-1; i++){
		for(j=i+1; j < NUM_OBJ; j++){
			//2d
			double dx = objList[i].pos[0] - objList[j].pos[0];
			double dy = objList[i].pos[1] - objList[j].pos[1];
			double dist = sqrt((dx*dx)+(dy*dy));

			if (dist < rad){
				double* ret = collide(objList[i].pos, objList[j].pos);
				print_collision(i, j, s);
				objList[i].pos[2] = ret[0];
				objList[i].pos[3] = ret[1];
				objList[j].pos[2] = ret[2];
				objList[j].pos[3] = ret[3];
			}
		}
				
	}
	print_step(objList, s);	
	return objList;
}

void print_collision(int a, int b, int d){
	int i;
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	fprintf(f,"Collision between object %d and %d at timestep %d.\n", a, b, d);
	fclose(f);
}

void print_timestep(int d) {
	int i;
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	fprintf(f,"--Time step %d--\n", d);	
	fclose(f);
}


double* collide(double* a, double* b){
	//obj a = *objList[i];
	//obj b = *objList[j];
	double* out = malloc(sizeof(double)*4);
	
	out[0] = b[2];
	out[1] = b[3];
	out[2] = a[2];
	out[3] = a[3];
	
	/*double enA = a[2]*a[3];
	double enB = b[2]*b[3];
	
	if (enA > enB){
		out[2] = (0 - (a[2] + b[2]))/2;
		out[3] = (0 - (a[3] + b[3]))/2;
		out[0] = (0 + (a[2] + b[2]))/2;
		out[1] = (0 + (a[3] + b[3]))/2;
	

		
	} else {
		out[0] = (0 - (a[2] + b[2]))/2;
		out[1] = (0 - (a[3] + b[3]))/2;
		out[2] = (0 + (a[3] + b[3]))/2;
		out[3] = (0 + (a[3] + b[3]))/2;
		
	}
	*/
	return out;
	//*objList[i] = a;
	//*objList[j] = b;*/
}

void print_step(obj* objList, int d){
	int i;
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	for(i=0; i < NUM_OBJ; i++){
		fprintf(f,"Object no: %d -{%.2f,%.2f,%.2f,%.2f} \n", i, objList[i].pos[0], objList[i].pos[1],
				objList[i].pos[2], objList[i].pos[3]);
	}
	fclose(f);
	
}

obj* gen_first_step(){
	obj* objList = malloc(sizeof(obj)*NUM_OBJ);
	int i;
	obj a;
	double* tA = malloc(sizeof(double)*4);
	//The following values are the external object / bullet / white ball
	tA[0] = -5;
	tA[1] = -2;
	tA[2] = 6;
	tA[3] = 5;
	a.pos = tA;
	objList[0] = a;
	for (i = 1; i < NUM_OBJ; i++){
		//The following values are the asteroid cluster
		obj t;
		double* temp = malloc(sizeof(double)*4);
		double r = (double) rand_int(X_MAX);
		double r2 = (double) rand_int(Y_MAX);
		double r3 = (double) (-rand_int(20)%X_VEC_MAX);
		double r4 = (double) (-rand_int(30)%Y_VEC_MAX);
		temp[0] = r;
		temp[1] = r2;
		temp[2] = r3;
		temp[3] = r4;
		t.pos = temp;
		objList[i] = t;
	}
	/*printf("First.a: %d,%d,%d,%d First.b: %d,%d,%d,%d\n", objList[0].pos[0], objList[0].pos[1],
				objList[0].pos[2], objList[0].pos[3],objList[1].pos[0], objList[1].pos[1], objList[1].pos[2], objList[1].pos[3]);
*/
	return objList;
}

/* Returns an integer in the range [0, n).
 *
 * Uses rand(), and so is affected-by/affects the same seed.
 */
int rand_int(int n) {
	//srand(time(NULL));
	int r = rand();
    return r % n;
}


double cpu_time ( void )

/*******************************************************************************/
/*
  Purpose:
     CPU_TIME reports the total CPU time for a program.

  Modified:
    27 September 2005

  Author:
    John Burkardt

  Parameters:
    Output, double CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}