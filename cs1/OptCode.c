/*
COSC3500 Assignment 1 - 2015
AMSED Validation Program
*/
#define DP

#ifdef SP
#define REAL float
#define ZERO 0.0
#define ONE 1.0
#define PREC "Single "
#endif

#ifdef DP
#define REAL double
#define ZERO 0.0e0
#define ONE 1.0e0
#define PREC "Double "
#endif

#define MATSIZE 700
#define VALIDATEFLG 1 //Set to zero and change MATSIZE to generate smaller random test matrices
#define FILENAME "TestData.c35"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

/*******************************************************************************/
/* Function Prototypes 														   */
/*******************************************************************************/
//**OPTIMISE THESE FUNCTIONS**
int CompRMat ( REAL R[][MATSIZE], REAL Vec1[], REAL Vec2[] );
int CompareVecMax ( REAL ResultMat[3][2], int n, REAL vec1[], REAL vec2[] );
int ExtractVariance( REAL PMat[][MATSIZE], REAL *VarVec);
void MultMatrix(REAL res[MATSIZE][MATSIZE], REAL a[MATSIZE][MATSIZE], REAL b[MATSIZE][MATSIZE]);
REAL madd(REAL a, REAL b, REAL c);
REAL ddot ( int n, REAL *dx, REAL *dy);
int LUDecomp( REAL* U, REAL* L, REAL* P, int n);
int Eye( REAL* A, int n); //1D version of Eye 2D
int SwapMatRows(REAL* A, const int R1, const int R2, const int n);
int Inv( REAL* A, int n);
int convert1DTo2DArray(REAL *sMat, REAL OutMat[MATSIZE][MATSIZE]);
int convert2DTo1DArray(REAL OutMat[][MATSIZE], REAL *sMat);
int Eye2D( REAL OutMat[][MATSIZE]);
int TransposeMatrix( REAL OutMat[][MATSIZE], REAL A[][MATSIZE]);
int addSubMatrix( REAL A[][MATSIZE], int Op, REAL B[][MATSIZE]);

//**TOOLS**
//	- Functions for you to use.
//	- Don't waste your time optimising these, but feel free to use them!
double cpu_time(void);
int disp( char* stext);
int PrintVec( REAL* vec, int n);
int PrintMatrix( REAL mat[MATSIZE][MATSIZE]);
REAL rand01();
void Pmatgen ( REAL P[][MATSIZE]);
void Hmatgen ( REAL H[][MATSIZE]);
void MeasErrorVecGen ( REAL *V1, REAL *V2);

/*******************************************************************************/

int main ( void )

/*******************************************************************************/
/*
  Purpose:
	1. Compute AMSED.
		- Adapted from a routine to compute the nth dimensional covariance time 
		  update of a Kalman filter.
	2. Validation of AMSED method.
*/
{
	REAL start, time;
	int val = MATSIZE;
	int *check_MatSize = &val;
	static REAL ResultMat[3][2], Pnew [MATSIZE][MATSIZE], I [MATSIZE][MATSIZE], K [MATSIZE][MATSIZE], P [MATSIZE][MATSIZE], H [MATSIZE][MATSIZE], R [MATSIZE][MATSIZE], TempMat1 [MATSIZE][MATSIZE], TempMat2 [MATSIZE][MATSIZE], S [MATSIZE][MATSIZE];
	static REAL S_A [MATSIZE*MATSIZE], V1[MATSIZE], V2[MATSIZE], PVar[MATSIZE], PnewVar[MATSIZE];
	FILE *fptr;

/*-----------------------------------------------------------------------------*/
/*      \/      Automatic Matrix Generation & Test Data File IO        \/      */
    if (VALIDATEFLG == 0){ //--------------Generate Random Matrices and Vectors
		// Generate P Matrix
		Pmatgen ( P );
		
		// Generate H Matrix
		Hmatgen ( H );
		
		// Generate Measurement Characteristic Vectors
		MeasErrorVecGen ( V1, V2);

	}else{ //----------------------------------Read Matrices & Vectors From File
		fptr = fopen(FILENAME, "rb");
		if (fptr!=NULL)
		{
			fread(check_MatSize, sizeof(int), 1, fptr);
			fread(P, sizeof(P), 1, fptr);
			fread(H, sizeof(H), 1, fptr);
			fread(V1, sizeof(V1), 1, fptr);
			fread(V2, sizeof(V2), 1, fptr);
			fclose(fptr);
		}else{
			printf("***Error! Failed to open:> %s\n", FILENAME);
			return;
		}
		
		//Check MATSIZE has been set correctly for validation mode
		if (MATSIZE != *check_MatSize){
				printf("***Error!***\nYou have attempted to read a file which was created using MATSIZE = %i\n", *check_MatSize);
				printf("Either:\n");
				printf("\t- change MATSIZE to %i and recompile before you rerun the program\n", *check_MatSize);
				printf("\tOR\n");
				printf("\t- change VALIDATEFLG to 0 and recompile before you rerun the program.\n");
				return 1;
		}
	}
/*      /\      Automatic Matrix Generation & Test Data File IO        /\      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	start = cpu_time();
	
	/*********************
	* Main Program Start *
	*********************/
	printf("Beginning Program: MatSize=%d, Precision=%s\n", MATSIZE, PREC);
	disp("Tasks:");
    
	//disp("V1:");
	//PrintVec(V1,MATSIZE);
	//disp("V2:");
	//PrintVec(V2,MATSIZE);
	
    //Create Measurement Error matrix R
    printf("\t- Computing R Matrix...");fflush(stdout);
		CompRMat ( R, V1, V2 );
	disp("complete.");
	
	//disp("R:");
	//PrintMatrix(R);
    
    //Evaluate S = H*P*H^T+R
	
	//----------- time measurement
	REAL startS = cpu_time();
	printf("\t- Computing S Matrix...");fflush(stdout);
		MultMatrix(TempMat1, H, P); 		//H*P
		TransposeMatrix( TempMat2, H); 		//H^T
		MultMatrix(S, TempMat1, TempMat2); 	//S = (K*P)*H^T
		addSubMatrix( S, 1, R); 			//S = (K*P*H^T)+R
	disp("complete.");
    
	REAL stopS = cpu_time();
	REAL difS = stopS-startS;
	printf("<Time of S: %.8f>\n", difS);
	
	//disp("S:");
	//PrintMatrix(S);
	REAL startK = cpu_time();
    //Evaluate K = P*H^T*S^-1
    printf("\t- Computing K Matrix...");fflush(stdout);
		TransposeMatrix( TempMat1, H); 		//H^T
		MultMatrix(TempMat2, P, TempMat1); 	//P*H^T
		convert2DTo1DArray(S,S_A);			//S 2d array incompatible with inv(), need to convert to SM
		Inv( S_A, MATSIZE);					//S^-1
		convert1DTo2DArray(S_A,S);			//Convert SM version of S^-1 back into 2D array
		MultMatrix(K, TempMat2, S);			//K=(P*H^T)*S^-1
	disp("complete.");
	
	REAL stopK = cpu_time();
	REAL difK = stopK-startK;
	printf("<Time of K: %.8f>\n", difK);
	//disp("K:");
	//PrintMatrix(K);
	
	REAL startP = cpu_time();
	//Evaluate Pnew = (I-K*H)P
	printf("\t- Computing Pnew Matrix...");fflush(stdout);
		Eye2D( I );							//I
		MultMatrix(TempMat1, K, H); 		//K*H
		addSubMatrix( I, -1,TempMat1); 		//I-K*H *NOTE* I = I-K*H
		MultMatrix(Pnew, I, P);				//(I-K*H)P
	disp("complete.");
	REAL stopP = cpu_time();
	REAL difP = stopP-startP;
	printf("<Time of P: %.8f>\n", difP);
    
	//disp("Pnew:");
	//PrintMatrix(Pnew);
	
	//Compare variance along original axes
	printf("\t- Comparing Variance...");fflush(stdout);
	ExtractVariance( P, PVar);
	ExtractVariance( Pnew, PnewVar);
	CompareVecMax ( ResultMat, MATSIZE, PVar, PnewVar );
	disp("complete.");
	
	//Report result to user
	printf("\t- Presenting Result.\n");
	if (ResultMat[0][0] == 1.0){
		printf("\tP contained the largest error along the reference axes:\n\t\tDimension: %i\n\t\tVariance: %f\n",(int)ResultMat[2][0]+1,ResultMat[1][0]);
		printf("\tPnew's largest error (along reference axes) was:\n\t\tDimension: %i\n\t\tVariance: %f\n",(int)ResultMat[2][1]+1,ResultMat[1][1]);
	}else{
		printf("\tPnew contained the largest error along the reference axes:\n\t\tDimension: %i\n\tVariance: %f\n",(int)ResultMat[2][1]+1,ResultMat[1][1]);
		printf("\tP's largest error (along reference axes) was:\n\t\tDimension: %i\n\t\tVariance: %f\n",(int)ResultMat[2][0]+1,ResultMat[1][0]);
	}

	/*******************
	* Main Program End *
	*******************/
	
	time = cpu_time() - start;
	printf("\n<Program completed in %.3f seconds.>\n\n", time);
	return 0;
}

/*******************************************************************************/

int CompRMat ( REAL R[][MATSIZE], REAL Vec1[], REAL Vec2[] )

/*******************************************************************************/
/*
  Purpose:
	Implements the measurement error matrix (R) initialisation function:
		R = diag(v1^2/dot(v1,v2)^2)
	
  Author:
    Tyler Hobson
*/
{
	int ii,row,col;
	REAL d12, d21, TempVec[MATSIZE];
	
	d12 = ddot( MATSIZE, Vec1, Vec2 );
	d21 = ddot( MATSIZE, Vec2, Vec1 );
	for (ii=0;ii < MATSIZE; ii++){
		 TempVec[ii] = Vec1[ii]*Vec1[ii] / (ddot( MATSIZE, Vec1, Vec2 )*ddot( MATSIZE, Vec2, Vec1 ));
	}
	
	for (row = 0;row<MATSIZE;row++){
		for (col = 0;col<MATSIZE;col++){
			if (row == col){
				R[row][col] = TempVec[row];
			}else{
				R[row][col] = ZERO;
			}
		}
	}
	
	return 0;

}
/*******************************************************************************/

int CompareVecMax ( REAL ResultMat[3][2], int n, REAL vec1[], REAL vec2[] )

/*******************************************************************************/
/*
  Purpose:
    CompareVecMax compares the maximum values contained in two vectors of length n.
	The results of the comparison are stored in a formatted 3x2 output matrix (ResultMat)
	as shown below:
	
	ResultMat=
		Vec1 Vec2
		 []  []  Larger Flag (value = 1 or 0) - indicates which vector is larger
		 []  []  Largest value from each vector
		 []  []  Dimension of largest value in each vector
	
  Author:
    Tyler Hobson
*/
{
	int ii;
	REAL maxValX = -99999999.9, maxValY = -99999999.9;

	//Find and record maximums of each vector
	for (ii=0;ii < n; ii++)
	{
		 if (vec1[ii] > maxValX){
			 maxValX = vec1[ii];
			 ResultMat[1][0] = maxValX;
		 }
	}
	for (ii=0;ii < n; ii++)
	{
		 if (vec2[ii] > maxValY){
			 maxValY = vec2[ii];
			 ResultMat[1][1] = maxValY;
		 }
	}

	//Find when the array with larger values is first larger than the other
	for (ii=0;ii < n; ii++)
	{
		if (maxValX > maxValY){
			if (ii == 0){
				ResultMat[0][0] = 1.0;
			}
			if (vec1[ii] == maxValX){
				ResultMat[2][0] = (REAL)ii;
			}
			if (vec2[ii] == maxValY){
				ResultMat[2][1] = (REAL)ii;
			}
		}else{
			if (ii == 0){
				ResultMat[0][1] = 1.0;
			}
			if (vec1[ii] == maxValX){
				ResultMat[2][0] = (REAL)ii;
			}
			if (vec2[ii] == maxValY){
				ResultMat[2][1] = (REAL)ii;
			}
		}
	}

	return 0;

}
/*******************************************************************************/

int ExtractVariance( REAL PMat[][MATSIZE], REAL *VarVec)

/*******************************************************************************/
/* 
  Purpose:
    Extracts the reference axis variances from a covariance matrix, 
	storing the result in the VarVec vector.
	
  Arguments:
		PMat 	- Input Matrix
		VarVec  - A vector containing the variances

  Author:
    Tyler Hobson

*/
{
	int row, col;
	REAL DSUM, MAX, SUM, NormVec[MATSIZE];
	
	//Prepare for matrix norm computation
	MAX = 0;
	for (col = 0; col < MATSIZE; col++){
		if (col > 0){
			DSUM = ddot(MATSIZE,(*PMat),((*(PMat)) + (MATSIZE * col)));
		}
		SUM = 0;
		for (row = 0; row < MATSIZE; row++){
			if (PMat[row][col] < 0){
				SUM = SUM + -PMat[row][col];
			}else{
				SUM = SUM + PMat[row][col];
			}
		}
		if (SUM > MAX){
			MAX = SUM;
		}
	}
	
	//Obtain main diagonal
	for (col = 0; col < MATSIZE; col++){
		for (row = 0; row < MATSIZE; row++){
			if (col==row){
				VarVec[row] = PMat[row][col];
			}
		}
	}
	
	return 0;
}
/*******************************************************************************/

void MultMatrix(REAL Res[][MATSIZE], REAL A[][MATSIZE], REAL B[][MATSIZE])

/*******************************************************************************/
{
/*
  Purpose:
    MultMatrix multiplies square matrices, A and B, storing the result in Res:
	Res = A*B
	
  Author:
    Tyler Hobson
*/
	REAL sum, temp, temp2, temp3, temp4, temp5;
	int i,j,k;
	
	for (i=0; i<MATSIZE; i++) {
		for (j=0; j<MATSIZE; j++) {
			sum = 0.0;
			for (k=0; k<MATSIZE;) {
				temp = A[i][k] * B[k][j];
				temp2 = A[i][k+1] * B[k+1][j];
				temp3 = A[i][k+2] * B[k+2][j];
				temp4 = A[i][k+3] * B[k+3][j];
				temp5 = A[i][k+4] * B[k+4][j];
				sum += temp + temp2 + temp3 + temp4 + temp5;
				k += 5;
				//sum = madd(sum,A[i][k],B[k][j]);
			}
			Res[i][j] = sum;
		}
	}
}
/*******************************************************************************/

REAL madd(REAL a, REAL b, REAL c)

/*******************************************************************************/
{
/*
  Purpose:
    When performing matrix and vector operations it is common to add and multiply
	values. Therefore to save time (writing code), madd will do this in one 
	command.
  Output:
	a + b * c
	
  Author:
    Tyler Hobson
*/
	return a + b * c;
}
/*******************************************************************************/

REAL ddot ( int n, REAL *dx, REAL *dy )

/*******************************************************************************/
/*
  Purpose:
    DDOT evaluates the dot product of two vectors dx & dy of length n.
	
  Output:
	Returns dot product.
	
  Author:
    Tyler Hobson
*/
{
	REAL dtemp, temp, SUM = 0.0;
	int i,ix,iy,m,mp1;
	
	//Scale Matrix by mp1
	dtemp = ZERO;
	
	//Perform dot product - Check for error condition
	if(n <= 0)
	{
		return(ZERO);
	}
	else if (n > 0)
	{
		for (i=0;i < n; i++)
		{
			temp = dx[i]*dy[i];
			SUM += temp;
			//SUM = madd(SUM,dx[i],dy[i]);
			if (i==n-1)
			{
				dtemp = SUM;
			}
		}
		return(dtemp);
	}
}
/*******************************************************************************/

int LUDecomp( REAL* U, REAL* L, REAL* P, int n)

/*******************************************************************************/
/* 
  Purpose:
    Determines the LU Decomposition of a square matrix A
    A gets replaced by U (i.e. set 1st argument = A; after execution A = U)

  Note:
	Function uses 1D array matrix notation
	
  Author:
    Tyler Hobson

*/
{
    
    int ErrFlag = 0;
    int rowMax;
    REAL l = 0.0;
    int rowcnt = n-1;
    int ii, col, row;
    
    
    Eye(L,n); //Ensure L = identity matrix
    for(ii = 0;ii<n*n;ii++) //Make P a copy of the identity matrix
        P[ii] = L[ii];
    

    //%Rearrange matrix for numeric stability
    for(col = 0;col<n-1;col++){
        //find row the maximum resides in
        
        rowMax = col;
        for (ii = col+1;ii<n;ii++){  //Find max
            if(U[n*col+ii] > U[n*col+rowMax]){
                rowMax = ii;
            }
        }

        if(rowMax != col){
            //want to swap row = col with rowMax
            SwapMatRows( U, col, rowMax, n);
            SwapMatRows( P, col, rowMax, n);
        }        
        
    }

    //%LU Decomposition
    col = 0;
    rowcnt = n-1;
    while (rowcnt > 0){
        if (U[col+col*n] == 0.0){
            for(ii = 0;ii<n*n;ii++) //Matrix is of low rank
            U[ii] = 0.0/0.0;
            ErrFlag = 1;
            return ErrFlag;
        }else{
            for (row = n-rowcnt;row<n;row++){
                    l = U[row+col*n]/U[col+col*n];
                    for (ii = 0;ii < n;ii++)
                        U[row+ii*n] = U[row+ii*n]-U[col+ii*n]*l;
                    L[row+col*n] = l;
            }
            rowcnt--;
            col++;
        }
    }
                
    return ErrFlag;
}
/*******************************************************************************/

int Eye( REAL* A, int n)

/*******************************************************************************/
/* 
  Purpose:
    Creates an nxn identity matrix

  Note:
	Function uses 1D array matrix notation
	
  Author:
    Tyler Hobson

*/
{
	int col, row;
    for (col = 0;col<n;col++){
        for (row = 0;row<n;row++){
            if(row == col){
               A[col*n+row] = 1.0;
            }else{
               A[col*n+row] = 0.0;
            }
        }
    }
	
	return 0;
}
/*******************************************************************************/

int Eye2D( REAL OutMat[][MATSIZE])

/*******************************************************************************/
/* 
  Purpose:
    Creates an nxn identity matrix

  Author:
    Tyler Hobson

*/
{
	int col, row;
    for (row = 0;row<MATSIZE;row++){
		for (col = 0;col<MATSIZE;col++){
            if(row == col){
               OutMat[row][col] = 1.0;
            }else{
               OutMat[row][col] = 0.0;
            }
        }
    }
	
	return 0;
}
/*******************************************************************************/

int SwapMatRows(REAL* A, const int R1, const int R2, const int n)
	
/*******************************************************************************/
/* 
  Purpose:
    Swaps 2 rows of a square matrix

  Note:
	Function uses 1D array matrix notation
	
  Author:
    Tyler Hobson

*/
{
    REAL TempVar;
    int col;
    
    for(col = 0;col<n;col++){       //for all columns
        TempVar = A[R1+col*n];      //Store R1 value
        A[R1+col*n] = A[R2+col*n];  //replace R1 with R2 value
        A[R2+col*n] = TempVar;      //replace R2 with stored value
    }

	return 0;
}
/*******************************************************************************/

int Inv( REAL* A, int n)

/*******************************************************************************/
/* 
  Purpose:
    Determines the Matrix inverse of an nxn matrix A, stores the result back into A
    Note: - Method does not handle low rank or ill condition matrices

  Note:
	Function uses 1D array matrix notation
	
  Author:
    Tyler Hobson

*/
{
    int ii,ErrFlag = 0;
    
//** First principals method - for small matrices   
    if (n < 4){
        REAL a [9];
        REAL Det;

        //Copy A into a
        for (ii = 0;ii<n*n;ii++){
            a[ii] = A[ii];
        }

        //method
        if (n == 2){
            Det = a[0]*a[3]-a[2]*a[1];
            A[0] = (1/Det)*a[3];
            A[1] = (1/Det)*-a[1];
            A[2] = (1/Det)*-a[2];
            A[3] = (1/Det)*a[0];
        }else if (n == 3){
            Det = a[0]*(a[8]*a[4]-a[5]*a[7])-a[1]*(a[8]*a[3]-a[5]*a[6])+a[2]*(a[7]*a[3]-a[4]*a[6]);
            A[0] = (1/Det)*(a[8]*a[4]-a[5]*a[7]);
            A[1] = (1/Det)*(-(a[8]*a[1]-a[2]*a[7]));
            A[2] = (1/Det)*(a[5]*a[1]-a[2]*a[4]);
            A[3] = (1/Det)*(-(a[8]*a[3]-a[5]*a[6]));
            A[4] = (1/Det)*(a[8]*a[0]-a[2]*a[6]);
            A[5] = (1/Det)*(-(a[5]*a[0]-a[2]*a[3]));
            A[6] = (1/Det)*(a[7]*a[3]-a[4]*a[6]);
            A[7] = (1/Det)*(-(a[7]*a[0]-a[1]*a[6]));
            A[8] = (1/Det)*(a[4]*a[0]-a[1]*a[3]);
        }
    }else{
//** LU Decomposition Method - Believed Numerically Unstable for near singular matrices
        
        int row,col;
        REAL DotResult;
		REAL *U, *L, *P, *Z;
		
		//Allocate memory for LU Decomp Matrices
		U = (REAL*) malloc (n*n*sizeof(REAL));
		L = (REAL*) malloc (n*n*sizeof(REAL));
		P = (REAL*) malloc (n*n*sizeof(REAL));
		Z = (REAL*) malloc (n*sizeof(REAL));

        if (n > 10000){
            ErrFlag = 1;
			free (U);
			free (L);
			free (P);
			free (Z);
            return ErrFlag;
        }
        for (ii = 0;ii<n*n;ii++)
            U[ii] = A[ii]; //Make a copy of A

        ErrFlag = LUDecomp(U,L,P,n);
        if (ErrFlag == 0){
            for (ii = 0;ii<n*n;ii++){
                A[ii] = 0.0; //Turn A into a matrix of zeros
            }

            //%ToDo: write an equivalent function to cond as a result integrity check

            for (col = 0;col<n;col++){
                for (ii=0;ii<n;ii++){ //%set/reset Z
                    Z[ii] = 0.0;
                }
                //%LZ = C : Using forward substitution
                for (row=0;row<n;row++){
                        DotResult = 0.0;
                        for (ii=0;ii<n;ii++)
                            DotResult = DotResult+L[row+ii*n]*Z[ii]; 
                        if (row == col){ //Hard code the identity matrix in to matrix B
                            Z[row] = (1.0-DotResult)/L[col+col*n];
                        }else{
                            Z[row] = -DotResult/L[col+col*n];
                        }   
                }

                //%UX = Z : Using backward substitution
                for (row = n-1;row >= 0;row--){
                    DotResult = 0.0;
                    for (ii=0;ii<n;ii++)
                        DotResult = DotResult+U[row+ii*n]*A[ii+col*n]; 
                    A[row+n*col] = (Z[row]-DotResult)/U[row+row*n];
                }
            }

            //%Correct for permutation - X = X*P (A = U*P)
            for (ii = 0;ii<n*n;ii++)
                U[ii] = A[ii]; //Make a copy of A
            MatMult(A,U,P,n);
        }
		
		free (U);
		free (L);
		free (P);
		free (Z);
    }
    return ErrFlag;
}
/*******************************************************************************/

int MatMult( REAL* C, REAL* A, REAL* B, int n)

/*******************************************************************************/
/* 
  Purpose:
    Determines the matrix product of two nxn matrices A and B, 
	storing the result in matrix C.

  Note:
	Function uses 1D array matrix notation
	
  Author:
    Tyler Hobson

*/
{
    REAL SUM;
	int row, col, ii;
    
    for (row=0;row<n;row++){
        for (col=0;col<n;col++){
            SUM = 0.0;
            for (ii=0;ii<n;ii++){
                SUM = SUM + A[row+ii*n]*B[col*n+ii];
            }
            C[row+col*n] = SUM;
        }
        
    }
	
	return 0;
}
/*******************************************************************************/

int convert1DTo2DArray( REAL *Mat1D, REAL Mat2D[][MATSIZE]){

/*******************************************************************************/
/* 
  Purpose:
    convert1DTo2DArray and convert2DTo1DArray allow you to change a matrix between
	a 2D and 1D array data type.
	
  Arguments:
		Mat1D 	   - Input 1D Matrix
		Mat2D 	   - Output 2D Matrix

  Author:
    Tyler Hobson

*/
	int row, col;
	
	for (col = 0; col < MATSIZE; col++){
		for (row = 0; row < MATSIZE; row++){
			Mat2D[row][col] = Mat1D[row*MATSIZE+col];
		}
	}
	
	return 0;
}
/*******************************************************************************/

int convert2DTo1DArray(REAL Mat2D[][MATSIZE], REAL *Mat1D){

/*******************************************************************************/
/* 
  Purpose:
    convert1DTo2DArray and convert2DTo1DArray allow you to change a matrix between
	a 2D and 1D array data type.
	
  Arguments:
		Mat2D 	   - Input 2D Matrix
		Mat1D 	   - Output 1D Matrix

  Author:
    Tyler Hobson

*/
	int row, col;
	
	for (row = 0; row < MATSIZE; row++){
		for (col = 0; col < MATSIZE; col++){
			Mat1D[row*MATSIZE+col] = Mat2D[row][col];
		}
	}
	
	return 0;
}
/*******************************************************************************/

int addSubMatrix( REAL A[][MATSIZE], int Op, REAL B[][MATSIZE]){

/*******************************************************************************/
/* 
  Purpose:
    Determines the sum or difference two nxn matrices A and B, 
	storing the result back into matrix A.
	
  Arguments:
		A 	   - Input Matrix
		B 	   - Input Matrix
		Op 	   - An integer indicating if A and B are to be added (Op >= 0) or subtracted (Op < 0)
		A 	   - The result of A<Op>B

  Author:
    Tyler Hobson

*/
	int row, col;
	
	for (row = 0; row < MATSIZE; row++){
		for (col = 0; col < MATSIZE; col++){
			if (Op < 0){
				A[row][col] = A[row][col]-B[row][col];
			}else{
				A[row][col] = A[row][col]+B[row][col];
			}
		}
	}
	
	return 0;
}
/*******************************************************************************/

int TransposeMatrix( REAL OutMat[][MATSIZE], REAL A[][MATSIZE]){

/*******************************************************************************/
/* 
  Purpose:
    Determines the transpose of an nxn matrix A, 
	storing the result in matrix OutMat.
	
  Arguments:
		A 	   - Input Matrix (matrix to transpose)
		OutMat - The result Transpose(A)

  Author:
    Tyler Hobson

*/
	int row, col;
	
	for (col = 0; col < MATSIZE; col++){
		for (row = 0; row < MATSIZE; row++){
			OutMat[row][col] = A[col][row];
		}
	}
	
	return 0;
}




/********************************************************************************
*===============================================================================*
*			    DO NOT OPTIMISE THE FUNCTIONS BELOW THIS BOX!!!!				*
*				     (But you are welcome to use them!!)						*
*===============================================================================*
********************************************************************************/
/*******************************************************************************/

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
/*******************************************************************************/

int disp( char* stext)

/*******************************************************************************/
/* 
  Purpose:
    Displays the text entered into this function to the screen and adds a newline char.

*/
{
     printf(stext);
	 printf("\n"); 
    
    return 0;
}
/*******************************************************************************/

int PrintVec( REAL *vec, int n)

/*******************************************************************************/
/* 
  Purpose:
    Prints the vector 'vec' to the screen.
	
  Inputs:
    vec - Array pointer
	n 	- Length of vector

*/
{
    int ii;
    
    for (ii=0;ii<n;ii++){
     printf("%f\n",vec[ii]);   
    }
    printf("\n");
    
    return 0;
}
/*******************************************************************************/

int PrintMatrix( REAL mat[][MATSIZE])

/*******************************************************************************/
/* 
  Purpose:
    Prints the matrix 'mat' to the screen

*/
{
    int ii,jj;
    
    for (ii=0;ii<MATSIZE;ii++){
		for (jj=0;jj<MATSIZE;jj++){
			printf("%0.2f ",mat[ii][jj]);   
		}
		printf("\n");
    }
	printf("\n"); 
    
    return 0;
}
/*******************************************************************************/

REAL rand01()

/*******************************************************************************/
/* 
  Purpose:
    Generates a number (REAL) between 0 and 1

*/
{
    return (REAL)rand() / (REAL)RAND_MAX ;
}
/*******************************************************************************/
//								DO NOT MODIFY!!								   //
void Pmatgen ( REAL P[][MATSIZE])
//								DO NOT MODIFY!!								   //
/*******************************************************************************/
/* 
  Purpose:
    Generates a valid Covariance Matrix (P) for testing

*/
{
  int ii, jj;
  
  srand((unsigned int)time(NULL)); //Set new random seed
  
  //Generate positive lower triangular matrix
  for ( ii = 0; ii < MATSIZE; ii++ )
  {
	  P[ii][ii] = (REAL)(rand() % 1000 + 1)+rand01(); //Main diagonal elements 1.99<=x<=1001
    for ( jj = ii+1; jj < MATSIZE; jj++ )
    {
        P[jj][ii] = rand01()*10.0; //Off-diagonal elements = pos(Norm(0 10))
    }
  }

  //Reflect lower triangular elements in upper triangular elements
  for ( jj = 1; jj < MATSIZE; jj++ )
  { 
    for ( ii = jj; ii < MATSIZE; ii++ )
    {
        P[jj-1][ii] = P[ii][jj-1];
    }
  }
  return;
}
/*******************************************************************************/
//								DO NOT MODIFY!!								   //
void Hmatgen ( REAL H[][MATSIZE])
//								DO NOT MODIFY!!								   //
/*******************************************************************************/
/* 
  Purpose:
    Generates a valid Observation Model Matrix (H) for testing

*/
{
  int ii, jj, row, col;
  
  srand((unsigned int)time(NULL)); //Set new random seed
  
  //Set elements to ZERO
  for ( row = 0; row < MATSIZE; row++ )
  {
    for ( col = 0; col < MATSIZE; col++ )
    {
        H[col][row] = ZERO;
    }
  }

  //Generate random direct associations between the state and measurement model
  for ( ii = 0; ii < MATSIZE; ii++ )
  { 
    if (ii%2 == 0){ //If ii is odd
		H[ii][ii] = 1;
	}
	for ( jj = ii+1; jj < MATSIZE; jj++ )
    {
		if ((rand() % 4 + 1) == 1){ // 1/4 chance 
			H[ii][jj] = 1;
		}
    }
  }
  return;
}
/*******************************************************************************/
//								DO NOT MODIFY!!								   //
void MeasErrorVecGen ( REAL *V1, REAL *V2)
//								DO NOT MODIFY!!								   //
/*******************************************************************************/
/* 
  Purpose:
    Generates measurement-error vectors for creating the R matrix

*/
{
  int ii;
  srand((unsigned int)time(NULL)); //Set new random seed
    
  for (ii=0;ii<MATSIZE;ii++){
	V1[ii] = rand01();
	V2[ii] = rand01()*0.01;
  }  
  
  return;
}