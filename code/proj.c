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
#include <mpi.h>

//Radius of objects
static double rad = 3;
//Filename to output to.
static char* filename = "output.txt";
const MPI_Comm comm=MPI_COMM_WORLD;
int myrank, size, collisions;
MPI_Datatype mpi_obj;

//Number of objects and timesteps to calculate
#define NUM_OBJ 10000
#define NUM_STEPS 5000
//X and Y max values for cloud objects
//#define X_MAX 50000
//#define Y_MAX 50000
//Max vector values for cloud objects (will be between 0 and max)
#define X_VEC_MAX 4
#define Y_VEC_MAX 4
//White ball initial state
#define W_XPOS -100
#define W_YPOS -100
#define W_XVEL 6
#define W_YVEL 6

/**An obj contains a array of floats of size 4
   index
   0 - x position
   1 - y position
   2 - x velocity
   3 - y velocity
*/
typedef struct {
	float xpos;
	float ypos;
	float xvel;
	float yvel;
} obj;

static void time_step(obj** stepCollec);
static obj* handle_step(obj* in, int s);
obj* gen_first_step();
int rand_int(int n, int proc);
double cpu_time(void);
void print_step(obj* step, int d);
void collide(int i, int j, obj* objList, obj* localObjList, int localNum);
void updateList(int start, int size, obj* objList, obj* localObjList);
void copyList(int start, int size, obj* objList, obj* localObjList);
void bufferCollide(int i, int j, obj* objList);
void updateBuffer(int* bufferBig, int* bufferNew, int size);
void updateBufferA(int* bufferA, int* bufferLocal, int size, int pos);
void print_init();
void compare(obj* objList, obj* localObjList);
void print_collision(int a, int b, int d);
void print_timestep(int d);
int X_MAX;
int Y_MAX;

void main(int argc, char** argv) {

	//Initialise MPI values
	MPI_Request Rrequests[2], Srequests[2];
    MPI_Status status[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &size);
	
	X_MAX = NUM_OBJ;
	Y_MAX = NUM_OBJ;
	
	if (myrank == 0) {
		printf("Started up...\n");
		//double start = cpu_time();
		remove(filename);
		print_init();
		//collisions = 0;
		//MPI_Bcast(&collisions, 1, MPI_INT, 0, comm);
    }
	collisions = 0;
	
	//Create MPI struct
	
	const int nitems=4;
    int          blocklengths[4] = {1,1,1,1};
    MPI_Datatype types[4] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    

    MPI_Aint offsets[4];

	offsets[0] = offsetof(obj, xpos);
	offsets[1] = offsetof(obj, ypos);
	offsets[2] = offsetof(obj, xvel);
	offsets[3] = offsetof(obj, yvel);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_obj);
    MPI_Type_commit(&mpi_obj);
	
	//Generate the objects in the initial time step
	obj* firstStep = malloc(sizeof(obj)*NUM_OBJ);
	if (myrank == 0){
		time_t seed;
		srand((unsigned) time(&seed));
		firstStep = gen_first_step();
		int t;
		for(t=1; t<size; t++){
			MPI_Send(&firstStep[0], NUM_OBJ, mpi_obj, t,5, comm);
		}
	} else {
		MPI_Status status;
		MPI_Recv(&firstStep[0], NUM_OBJ, mpi_obj, 0, 5, comm, &status);
	}
	MPI_Barrier;
		

	
	//Allocate memory for the collection of time steps
	obj** stepCollec = malloc(sizeof(obj)*NUM_OBJ);
	stepCollec[0] = firstStep;

	//Perform time stepping
	time_step(stepCollec);
	
	if (myrank != 0 ){
		MPI_Send(&collisions, 1, MPI_INT, 0, 4, comm);
	}

	
	if (myrank == 0) {
		int t, collisionsT;
		MPI_Status status[size];
		for(t=1; t<size; t++){
			MPI_Recv(&collisionsT, 1, MPI_INT, t, 4, comm, &status[t]);
			//printf("CollisionsT %d \n", collisionsT);
			collisions += collisionsT;
		}

		printf("Collisions %d \n", collisions);
	}

	
	//double stop = cpu_time()-start;
	
	//printf("Finished in %f seconds. Exiting...\n", stop);
	MPI_Finalize();
	return;
}
/*Generate initial state. Object a is the whiteball/bullet
* Objects 1 - NUM_OBJ are cloud objects defined in top of file
*/
obj* gen_first_step(){
	obj* objList = malloc(sizeof(obj)*NUM_OBJ);
	int i;
	obj a;
	//The following values are the external object / bullet / white ball
	a.xpos = W_XPOS;
	a.ypos = W_YPOS;
	a.xvel = W_XVEL;
	a.yvel = W_YVEL;
	objList[0] = a;
	#pragma omp parallel shared(objList)
	{
		int proc = omp_get_thread_num();
		#pragma omp for 
		for (i = 1; i < NUM_OBJ; i++){
			//The following values are the asteroid cluster
			obj t;
			//Randomly generate each object (object conists: xpos, ypos, xve, yvel)
			float r = (float) rand_int(X_MAX, proc);
			float r2 = (float) rand_int(Y_MAX , proc);
			float r3 = (float) (-rand_int(20, proc)%X_VEC_MAX);
			float r4 = (float) (-rand_int(30, proc)%Y_VEC_MAX);
			t.xpos = r;
			t.ypos = r2;
			t.xvel = r3;
			t.yvel = r4;
			objList[i] = t;
		}
	}
	return objList;
}

/*Perform time steps for whole collection
*/
void time_step(obj** stepCollec) {
	int i;
	for(i=1; i < NUM_STEPS+1; i++){
		obj* step = handle_step(stepCollec[0], i);
		//free(stepCollec);
		stepCollec[0] = step;
		//free(stepCollec[i-1]);
	}
	
	return;
	
}
/*Perform individual time steps
*/
obj* handle_step(obj* in, int s) {
	int i,j;
	if (myrank == 0){
		print_timestep(s);
	}
	
	obj* objList = malloc(sizeof(obj)*2*NUM_OBJ);
	#pragma omp parallel shared(objList)
	{
		//Do time advancement on all objects in timestep
		#pragma omp for 
		for(i=0; i < NUM_OBJ; i++){
			float tAx = in[i].xpos + in[i].xvel;
			float tAy = in[i].ypos + in[i].yvel;
			obj a;
			a.xpos = tAx;
			a.ypos = tAy;
			a.xvel = in[i].xvel;
			a.yvel = in[i].yvel;
			objList[i] = a;

		}
	}
	
	//Initalise thread unique domain decomp. values
	int localSize = NUM_OBJ/size;
	int localNum = myrank*localSize;
	int localTopVal = localNum + localSize;
	
	//Initalise buffer which stores collisions outside of thread domain
	int *buffer = malloc(sizeof(int)*2*NUM_OBJ*localSize);
	int bufferUsed = 0;
	
	//Initialise local (this thread) object list
	obj* localObjList = malloc(sizeof(obj)*2*localSize);
	copyList(localNum, localSize, objList, localObjList);
	
	int num = NUM_OBJ-1;
	//Check collisions between every object
	for(i=0; i < num; i++){
		for(j=localNum+1; j < localTopVal; j++){
			//If not the same object.
			if (i != j){
				
				float dx = objList[i].xpos - objList[j].xpos;
				float dy = objList[i].ypos - objList[j].ypos;
				float dist = sqrt((dx*dx)+(dy*dy));

				//If collision detected between objects (radius defined in header)
				if (dist < rad){
					

					//If not in thread domain, add object numbers to buffer for later
					if (i < localNum || i >= localTopVal){
						buffer[bufferUsed] = i;
						buffer[bufferUsed+1] = j;
						bufferUsed += 2;
					//Else do collision
					} else{
						collisions += 1;
						collide(i, j, objList, localObjList, localNum);
					}
				}
			}
		}
	}
	
	//Resize buffer to actual bufferUsed size
	int *newBuffer = malloc(sizeof(int)*2*bufferUsed);
	updateBuffer(buffer, newBuffer, bufferUsed*2);

	//Initalise MPI stuffs
	int nsends=0;
	MPI_Status status[size*4];
	MPI_Request Srequests[4], Rrequests[size*4];
	//If not root send local object list, buffer size and buffer to root
	if (myrank != 0 ){
		int bufferSize = bufferUsed;
		MPI_Send(&localObjList[0], localSize, mpi_obj, 0, 1, comm);
		nsends++;
		MPI_Send(&bufferSize, 1, MPI_INT, 0, 2, comm);
		nsends++;
		MPI_Send(&newBuffer[0], bufferSize, MPI_INT, 0, 3, comm);
		nsends++;

	//If root thread, recieve all updates
	} else {
	//Master thread receives all data
	//for each message
		// Update data on root
		// Update buffered collisions
		// Broadcast results for timestep
		
		
		//Initalise messy buffer things
		//BufferA is massive flat 2D buffer array, storing all the buffers
		//BufferUsedA is array storing all the buffer sizes
		//Buffer local is the current buffer being unloaded from MPI
		int bufferMult = 2*NUM_OBJ*localSize;
		int *bufferA = malloc(sizeof(int)*bufferMult*size);
		int *bufferUsedA = malloc(sizeof(int)*size);
		int* bufferLocal = malloc(sizeof(int)*bufferMult);

		//Update buffer from root thread
		updateBufferA(bufferA, buffer, bufferUsed, 0);
		int bufferPos = bufferUsed;
		bufferUsedA[0] = bufferUsed;

		//Update data from root thread
		updateList(0, localSize, objList, localObjList);
		
		int t;
		//For each thread, recieve data and update into buffer array and global object list
		for(t=1; t<size; t++){
			int tag, collisionsT;
			int localNum = localSize*t;
			int *localBuffer = malloc(sizeof(int)*2*NUM_OBJ);
			obj* localObjList = malloc(sizeof(obj)*localSize);
			int bufferUsed, b;
			MPI_Recv(&localObjList[0], localSize, mpi_obj, t, 1, comm, &status[t]);
			MPI_Recv(&bufferUsed, 1, MPI_INT, t, 2, comm, &status[t+1]);
			MPI_Recv(&bufferLocal[0], bufferUsed, MPI_INT, t, 3, comm, &status[t+2]);
			
			//Insert buffer array size
			bufferUsedA[t] = bufferUsed;
			
			//Update buffer array with buffer from thread
			updateBufferA(bufferA, bufferLocal, bufferUsed, bufferPos);
			
			//Increment position within flat 2D big buffer array
			bufferPos += bufferUsed;
			
			//Update global object list
			updateList(localNum, localSize, objList, localObjList);
			
			free(localBuffer);
			free(localObjList);
			
		}
		
		//Now all data has been loaded into root,
		//Update global object list with buffered collisions!
		bufferPos = 0;
		//For each thread/rank
		for(t=0; t<size; t++){
			int b;
			int bufferUsedLocal = bufferUsedA[t];
			//For each object pair
			for (b = 0; b < bufferUsedLocal; b++){
					collisions += 1;
					bufferCollide(bufferA[bufferPos+b], bufferA[bufferPos+b+1], objList);
					b += 1;
			}
			bufferPos += bufferUsedLocal;
		}
		//Broadcoast timestep results and print
		MPI_Bcast(&objList[0], NUM_OBJ, mpi_obj, 0, comm);
		print_step(objList, s);	
		free(bufferA);
		free(bufferUsedA);
		free(bufferLocal);
	}

	//Wait for all threads
	MPI_Barrier(comm);
	
	free(localObjList);
	free(newBuffer);
	free(buffer);

	return objList;
}

/*Method to compare local and objList (used for verifcation when using only 1 thread)
*/
void compare(obj* objList, obj* localObjList){
	remove("compare.txt");
	FILE *f = fopen("compare.txt", "a");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}	
	fclose(f);
	int i;
	for (i = 0; i < NUM_OBJ; i++){
		if (objList[i].xvel != localObjList[i].xvel || objList[i].yvel != localObjList[i].yvel){
			fprintf(f,"Diff found! %d %d - %d %d", objList[i].xpos, objList[i].ypos, localObjList[i].xpos, localObjList[i].ypos);
		}
	}	
}

/*Update used when resizing local buffer
*/
void updateBuffer(int* bufferBig, int* bufferNew, int size){
	int i;
	for (i = 0; i < size; i++){
		bufferNew[i] = bufferBig[i];
		bufferNew[i+1] = bufferBig[i+1];
		i++;
	}	
}

/*Update used when inserting thread buffers into buffer array
*/
void updateBufferA(int* bufferA, int* bufferLocal, int size, int pos){
	int i;
	for (i = 0; i < size; i++){
		bufferA[pos+i] = bufferLocal[i];
		bufferA[pos+i+1] = bufferLocal[i+1];
		i++;
	}	
}

/*Update object list with local list contents
*/
void updateList(int start, int size, obj* objList, obj* localObjList){
	int i;
	for (i = 0; i < size; i++){
		objList[start+i] = localObjList[i];
	}	
}

/*Copy contents of objList into local list
*/
void copyList(int start, int size, obj* objList, obj* localObjList){
	int i;
	for (i = 0; i < size; i++){
		localObjList[i] = objList[start+i];
	}	
}


/*Handle collision velocity calcualtions of two objects i and j.
* 
*/
void collide(int i, int j, obj* objList, obj* localObjList, int localNum){
	
	//int localTopVal = localNum + localSize;
	int iLocal = i-localNum;
	int jLocal = j-localNum;
	
	localObjList[iLocal].xvel = objList[j].xvel;
	localObjList[iLocal].yvel = objList[j].yvel;
	localObjList[jLocal].xvel = objList[i].xvel;
	localObjList[jLocal].yvel = objList[i].yvel;
}

/*Handle collision velocity calcualtions of two objects i & j, 
* called within buffer corrections
*/
void bufferCollide(int i, int j, obj* objList){
	
	//print_collision(i, j, 0);
	obj oi = objList[i];
	obj oj = objList[j];
	
	objList[i].xvel = oj.xvel;
	objList[i].yvel = oj.yvel;
	objList[j].xvel = oi.xvel;
	objList[j].yvel = oi.yvel;
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


/*Print timestep number to file.
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
		fprintf(f,"Object: %d - %.2f %.2f %.2f %.2f \n", i, objList[i].xpos, objList[i].ypos,
				objList[i].xvel, objList[i].yvel);
	}
	fclose(f);
	
}

/* Returns an integer in the range [0, n).
 */
int rand_int(int n, int proc) {
	//int r =  rand_r(&proc);
	int r =  rand();
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
