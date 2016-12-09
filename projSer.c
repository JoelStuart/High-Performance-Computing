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
*
* Program outputs to static filename defined below. Read this output
* for verification
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <mpi.h>

//Radius of objects
static double rad = 3;
//Filename to output to.
static char* filename = "outputSer.txt";
//const MPI_Comm comm=MPI_COMM_WORLD;
int myrank, size, collisions;

//Number of objects and timesteps to calculate
#define NUM_OBJ 5000
#define NUM_STEPS 50
//X and Y max values for cloud objects
#define X_MAX 5000
#define Y_MAX 5000
//Max vector values for cloud objects (will be between 0 and max)
#define X_VEC_MAX 4
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
double* collide(double* a, double* b);
void print_init();
void print_collision(int a, int b, int d);
void print_timestep(int d);

void main(void) {

	printf("Started up...\n");
	//double start = cpu_time();
	remove(filename);
	print_init();
	collisions = 0;
	//Generate the objects in the initial time step
	obj* firstStep = gen_first_step();
	
	//Allocate memory for the collection of time steps
	obj** stepCollec = malloc(sizeof(obj[NUM_OBJ])*NUM_STEPS);
	stepCollec[0] = firstStep;

	//Perform time stepping
	time_step(stepCollec);
	
	//double stop = cpu_time()-start;
	
	printf("Collisions: %d\n", collisions);


  return;
}
/*Generate initial state. Object a is the whiteball/bullet
* Objects 1 - NUM_OBJ are cloud objects defined in top of file
*/
obj* gen_first_step(){
	obj* objList = malloc(sizeof(obj)*NUM_OBJ);
	int i;
	obj a;
	double* tA = malloc(sizeof(double)*4);
	//The following values are the external object / bullet / white ball
	tA[0] = -100;
	tA[1] = -100;
	tA[2] = 6;
	tA[3] = 6;
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
	return objList;
}

/*Perform time steps for whole collection
*/
void time_step(obj** stepCollec) {
	int i;
	for(i=1; i < NUM_STEPS+1; i++){
		stepCollec[i] = handle_step(stepCollec[i-1], i);
		//free(stepCollec[i-1]);
	}
	
	return;
	
}
/*Perform individual time steps
*/
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
				//Ret is a holder array with x and y velocity values
				//For two objects i and j
				//Stored 0 - a.x, 1 - a.y, 2 - b.x, 3 - b.y
				//print_collision(i, j, s);
				collisions += 1;
				double* ret = collide(objList[i].pos, objList[j].pos);
				//print_collision(i, j, s);
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

/*Handle collision velocity calcualtions of two objects a and b.
* store results in return array out. {a.x, a.y, b.x, b.y}
*/
double* collide(double* a, double* b){
	double* out = malloc(sizeof(double)*4);
	//Reverse direction of travel vector
	out[0] = b[2];
	out[1] = b[3];
	out[2] = a[2];
	out[3] = a[3];
	return out;
}

/*Print collision occurance to file
*/
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

/*Print number of objects and timesteps to file
*/
void print_init() {
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}	
	fprintf(f,"%d %d\n", NUM_OBJ, NUM_STEPS);	
	fclose(f);
}


/*Print timestep number to file. Could have concievably written a general
* string print to file by now, but copy paste was quicker.
*/
void print_timestep(int d) {
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	//fprintf(f,"--Time step %d--\n", d);	
	fprintf(f,"TS %d\n", d);	
	fclose(f);
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
		//fprintf(f,"Object no: %d -{%.2f,%.2f,%.2f,%.2f} \n", i, objList[i].xpos, objList[i].ypos,
		//		objList[i].xvel, objList[i].yvel);
		fprintf(f,"Object: %d - %.2f %.2f %.2f %.2f \n", i, objList[i].pos[0], objList[i].pos[1],
				objList[i].pos[2], objList[i].pos[3]);
	}
	fclose(f);
	
}

/* Returns an integer in the range [0, n).
 */
int rand_int(int n) {
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