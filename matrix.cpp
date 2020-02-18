#include "matrix.hpp"
#include <cstdint>
#include <iostream>
 using namespace std;
/*
    Nous allons stocker des matrices dans des tablaux à 1 dimension:
    Par exemple, matrice de taille n x m est stockué comme:
        A_tab = [A_11, A_12, ... A_1m, A_21, ... A_2m, ... , A_1n, ..., A_nm]

    Par conséquant, un élement A_ij aurait quelle indice dans le tableau A_tab ? 

 */

/* Memory allocation for a matrix of size n x m and initilization to 0  */
double *allocateMatrix(uint64_t n,uint64_t m) {
  double *A;
  A = (double *) calloc (n * m, sizeof(double));
  return A;
}


/* Frees the memory allocated to matrix A
*/
void freeMatrix(double *A) {
    free(A);
}

/* Allocates a n sized vector and initializes all entries to 0 
*/
double *allocateVector(uint64_t n) {
  double *v; 
  v = (double *) calloc(n, sizeof(double));
  return v;
}

/* Trees the memory allocated to a vector
*/
void freeVector(double *v) {
  free(v);
}


/* Sets a n * m matrix A to all zeros */
void setMatrixZero(double *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
        /* Note that for a n x m matrix flattened to a 1D array, 
        element A_ij has index i * m + j
        */
      A[i * m + j] = 0.0; 
    }
  }
}

/* Sets a n * n matrix A to identity */
void setMatrixIdentity(double *A, uint64_t n) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
     A[i * n + j] = 0.0;
    }
    A[i * n + i] = 1.0;
  }
}



/* Copies a matrix  
*/
void copyMatrix(double *B, double *A, uint64_t n, uint64_t m) {
  uint64_t i,j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      B[i * m + j] = A[i * m + j]; 
    }
  }
}

/*
Writes a matrix to a stream. For example, writing a matrix to standard output is
writeMatrix(stdout, A, n, m);
A sream can also be a file. 
*/
void writeMatrix(FILE *stream, double *A, uint64_t n, uint64_t m)
{
	fprintf(stream, "%d %d \n", (int)n, (int)m);
	int i, j;
	for(i = 0; i < n; ++i)
	{
	      for(j = 0; j < m; ++j)
	      {
		      fprintf(stream, "%f \t", A[i * m + j]);
	      }
	      fprintf(stream, "\n");
	}
}



//The function computes the element-by-element abs of matrix A
void absMatrix(double *Aabs,double *A, uint64_t n, uint64_t m)
{
	uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            Aabs[i*m + j] = fabs(Aabs[i*m + j]);
		}
	}

}


/*
Performs addition of two matrix A (size n x m) and B (size n x m).
The result S = A + B is a n x m matrix.
We consider that S is allocated outside the function.
*/
void matrixAdd(double *S, double *A, double *B, uint64_t n, uint64_t m){
    uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            S[i*m + j] = A[i*m + j] + B[i*m + j];
		}
	}
}

/*
Performs subtraction of two matrix A (size n x m) and B (size n x m).
The result S = A - B is a n x m matrix.
We consider that S is allocated outside the function.
*/
void matrixSub(double *S, double *A, double *B, uint64_t n, uint64_t m){
    uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            S[i*m + j] = A[i*m + j] - B[i*m + j];
		}
	}
}



/* For a double m x n matrix A the function returns its maximum in absolute value
element. */
double getMaxInMatrix(double max, double *A, uint64_t n, uint64_t m)
{
	double maxA = fabs(A[0]);
	double current = fabs(A[0]);
	int i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
			current = fabs(A[i * m + j]);
			if(current > maxA)
				maxA = current;
		}
	}
 //element A_ij has index i * m + j

}


/* Rajouter les prototypes de vos méthodes ici. Par exemple */

/* Performs naive multiplication of matrix A (size p x k) by a matrix B (size k x r).
The result matrix S = A*B  is of size (p x r).
We assume that S has already been allocated outside the function.
*/
void        matrixMultiplyNaive (double *S, double *A, double *B, uint64_t p, uint64_t k, uint64_t r){
    

    for (int s = 0; s < r; s++){ 

      for (int i =0; i < p; i++){ 

        for (int j=0; j < k; j++){
            S[i*r+s] += A[i*k+j] * B[j*k+s];
        }

      }

    }
    
}

/* Performs a multiplication of two sqaure matrices A and B (size n x n) by Strassen algorithm.
    We assume that S has already been allocated outside the function.
*/
void        matrixMultiplyStrassen (double *S, double *A, double *B, uint64_t n){

    if(n == 2){
      double m1,m2,m3,m4,m5,m6,m7;
      m1 = (A[0] + A[3])*(B[0]+ B[3]);
      m2 = (A[2] + A[3])*B[0];
      m3 = (B[1] - B[3])*A[0];
      m4 = (B[2] - B[0])*A[3];
      m5 = (A[0] + A[1])*B[3];
      m6 = (A[2] - A[0])*(B[0]+ B[1]);
      m7 = (A[1] - A[3])*(B[2]+ B[3]);

      S[0] = m1+m4-m5+m7;
      S[1] = m3+m5;
      S[2] = m2 + m4;
      S[3] = m1 - m2 + m3 + m6;



    }else{
      double* tmp,*tmp2;
      tmp = allocateMatrix(n/2,n/2);
      tmp2 = allocateMatrix(n/2,n/2);

      double* resBlock0,*resBlock1,*resBlock2,*resBlock3;

      resBlock0 = allocateMatrix(n/2,n/2);
      resBlock1 = allocateMatrix(n/2,n/2);
      resBlock2 = allocateMatrix(n/2,n/2);
      resBlock3 = allocateMatrix(n/2,n/2);
      
      

      double* m1,*m2,*m3,*m4,*m5,*m6,*m7;
      m1 = allocateMatrix(n/2,n/2);
      m2 = allocateMatrix(n/2,n/2);
      m3 = allocateMatrix(n/2,n/2);
      m4 = allocateMatrix(n/2,n/2);
      m5 = allocateMatrix(n/2,n/2);
      m6 = allocateMatrix(n/2,n/2);
      m7 = allocateMatrix(n/2,n/2);

      
      double* A0 = allocateMatrix(n/2,n/2);
      double* A1 = allocateMatrix(n/2,n/2);
      double* A2 = allocateMatrix(n/2,n/2);
      double* A3 = allocateMatrix(n/2,n/2);
      double* B0 = allocateMatrix(n/2,n/2);
      double* B1 = allocateMatrix(n/2,n/2);
      double* B2 = allocateMatrix(n/2,n/2);
      double* B3 = allocateMatrix(n/2,n/2);

      //dividing the matrix into four smaller matrices
      for(int i = 0 ;i <= n-1; i++){
       if(i < n/2){
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             A0[(i * (n/2)) + j] = A[(i * n) + j ] ;
             B0[(i * (n/2)) + j] = B[(i * n) + j ] ;
           }else{
             A1[(i * (n/2)) + (j % (n/2))] = A[(i * n) + j ] ;
             B1[(i * (n/2)) + (j % (n/2))] = B[(i * n) + j ] ;
           }
         }
       }else{
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             A2[((i%2) * (n/2)) + j] = A[(i * n) + j ] ;
             B2[((i%2) * (n/2)) + j] = B[(i * n) + j ] ;
           }else{
             A3[((i%2) * (n/2)) + (j % (n/2))] = A[(i * n) + j ] ;
             B3[((i%2) * (n/2)) + (j % (n/2))] = B[(i * n) + j ] ;
           }
         }
       }
      

      }

     

      //m1
      matrixAdd(tmp , A0 , A3 , n/2 , n/2);
      matrixAdd(tmp2 , B0 , B3 , n/2 , n/2);
      matrixMultiplyStrassen(m1,tmp,tmp2,n/2);
      

      //m2
      matrixAdd(tmp , A2 , A3 , n/2 , n/2);
      matrixMultiplyStrassen(m2,tmp,B0,n/2);

      //m3
      matrixSub(tmp , B1 , B3 , n/2 , n/2);
      matrixMultiplyStrassen(m3,A0,tmp,n/2);

      //m4
      matrixSub(tmp , B2 , B0 , n/2 , n/2);
      matrixMultiplyStrassen(m4,A3,tmp,n/2);


      //m5
      matrixAdd(tmp , A0 , A1 , n/2 , n/2);
      matrixMultiplyStrassen(m5,tmp,B3,n/2);

      //m6
      matrixSub(tmp , A2 , A0 , n/2 , n/2);
      matrixAdd(tmp2 , B0 , B1 , n/2 , n/2);
      matrixMultiplyStrassen(m6,tmp,tmp2,n/2);

      //m7
      matrixSub(tmp , A1 , A3 , n/2 , n/2);
      matrixAdd(tmp2 , B2 , B3 , n/2 , n/2);
      matrixMultiplyStrassen(m7,tmp,tmp2,n/2);


      matrixAdd(tmp,m1,m4,n/2,n/2);
      matrixSub(tmp2,m7,m5,n/2,n/2);
      matrixAdd(resBlock0,tmp,tmp2,n/2,n/2);

      matrixAdd(resBlock1,m3,m5,n/2,n/2);

 
      matrixAdd(resBlock2,m2,m4,n/2,n/2);


      matrixAdd(tmp,m3,m6,n/2,n/2);
      matrixSub(tmp2,m1,m2,n/2,n/2);
      matrixAdd(resBlock3,tmp,tmp2,n/2,n/2);

      //filling the result matrix S
      for(int i = 0 ;i <= n-1; i++){
       if(i < n/2){
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             S[(i * n) + j ] =  resBlock0[(i * (n/2)) + j] ;
            
           }else{
             S[(i * n) + j ] = resBlock1[(i * (n/2)) + (j % (n/2))] ;
             
           }
         }
       }else{
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             S[(i * n) + j ] = resBlock2[((i%2) * (n/2)) + j]  ;
             
           }else{
             S[(i * n) + j ] = resBlock3[((i%2) * (n/2)) + (j % (n/2))];
             
           }
         }
       }
      

      }
      

    }

    
}

/* 
    Solves a system of linear equations Ax=b for a double-precision matrix A (size n x n).
    Uses iterative ascension algorithm. 
    After the procedure, x contains the solution of Ax=b.
    We assume that x has been allocated outside the function.
*/
void        SolveTriangularSystemUP   (double *x, double *A, double *b, uint64_t n){

    for(int i = n-1 ; i >= 0 ; i-- ){
      double pivot = A[(i * n) + i];
      b[i] = b[i] / pivot;
      x[i] = b[i];
      A[(i * n) + i] = pivot / pivot;

      for(int j = i-n ; j >= 0; j=j-n){
        A[(i * n) + j] = A[(i * n) + j] - (A[(i * n) + j]*pivot);
        b[j] = b[j] - (b[j]*b[i]);



      }
      



    }
}

/* 
    Performs Gauss elimination for given a matrix A (size n x n) and a vector b (size n).
    Modifies directly matrix A and vector b.
    In the end of the procedure, A is upper truangular and b is modified accordingly.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is impossible to triangularize. 
*/
bool        Triangularize           (double *A, double *b, uint64_t n){
  int toPivot = n +1;

  for(int i = 0 ; i <= n*n ; i = i + toPivot){

    if(A[i] == 0){  //swapping

      int j = i;

      while(A[j] == 0 && j < n*n){
        j = j + n;
      }

      if(j > n*n){
        return false;
      }

      for(int k = j%n ; k < n ; k++){  //swapping the lines
        double tmp = A[i + k] ;
        A[i+k] = A[j+k];
        A[j+k] = tmp;


      }

      int LinePivot = (int)(i / n);
      int indVect = (int)(j / n);  //index in the vector

      //swapping the lines in the vector
      double tmp = b[LinePivot] ;
      b[LinePivot] = b[indVect];
      b[indVect] = tmp;

    }
  

    //triangularizing
    for(int p = i+n; p <= n*n; p = p + n){

      double valUnderPivot = A[p];
      int lastInLine = n - (p%n);
      for(int y = 0; y < lastInLine ; y++){
        A[p+y] = A[p+y] - ((valUnderPivot/A[i])*A[i+y]);
      }

      int LinePivot = (int)(i / n);
      int indVect = (int)(p / n);  //index in the vector

      b[indVect] = b[indVect] - (valUnderPivot/A[i])*b[LinePivot];  //applying the same operation the the vector

    }



  }
    
    

    return true;
}

/*
    Solves a system of linear equations Ax=b, given a matrix A (size n x n) and vector b(size n).
    Uses Gauss elimination algorithm based on truangularization and the ascension solving.
    After the procedure, vector x contains the solution to Ax=b.
    We assume that x has been allocated outside the function.
        Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is of rank <n .
*/
bool        SolveSystemGauss        (double *x, double *A, double *b, uint64_t n){
    
    bool triangular = Triangularize(A,b,n);
    writeMatrix(stdout, A, n, n);
    if(triangular){
      SolveTriangularSystemUP(x,A,b,n);
      return true;
    }else{
      return false;
    }

    
}

