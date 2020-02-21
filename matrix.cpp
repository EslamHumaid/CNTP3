/* Authors: Eslam HUMAID, Abrahim BAMATRAF.
    Groupe : 485L
*/
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
    

    for (int s = 0; s < r; s++){ //the column of B.

      for (int i =0; i < p; i++){  //the line in A.

        for (int j=0; j < k; j++){  //the line in B and column in A.
            S[i*r+s] += A[i*k+j] * B[j*r+s]; //using the rule (element A_ij has index i * m + j).
        }

      }

    }
    
}

/*
determines wether a number can be written in the form 2^n
 */
bool isPowerOftwo(uint64_t n){
  if(n>0)
    {
        while(n%2 == 0)
        {
            n/=2;
        }
        if(n == 1)
        {
            return true;
        }
    }
    if(n == 0 || n != 1)
    {
        return false;
    }
}

/* Performs a multiplication of two sqaure matrices A and B (size n x n) by Strassen algorithm.
    We assume that S has already been allocated outside the function.
    
 */
void        matrixMultiplyStrassen (double *S, double *A, double *B, uint64_t n){
    if(n == 1){
      S[0] = A[0] * B[0] ;

    }else if(n == 2){ //in this case the matrix is already in the form 2x2 so no need to divide it into smaller ones.
      
      
      //using strassen formulas.
      double m1,m2,m3,m4,m5,m6,m7;
      m1 = (A[0] + A[3])*(B[0]+ B[3]);
      m2 = (A[2] + A[3])*B[0];
      m3 = (B[1] - B[3])*A[0];
      m4 = (B[2] - B[0])*A[3];
      m5 = (A[0] + A[1])*B[3];
      m6 = (A[2] - A[0])*(B[0]+ B[1]);
      m7 = (A[1] - A[3])*(B[2]+ B[3]);

      //the result matrix.
      S[0] = m1+m4-m5+m7;
      S[1] = m3+m5;
      S[2] = m2 + m4;
      S[3] = m1 - m2 + m3 + m6;



    }else{ // dividing the matrices into four (n/2)x(n/2) matrices.

      //Copies of A and B to change their dimensions if needed.
      double* copyA = A;  
      double* copyB = B;
      
      if(!isPowerOftwo(n)){  //if the dimension is not 2^n*2^n
        uint64_t validN = n; //the right dimension in the form 2^n.

        while(!isPowerOftwo(validN)){ //to find the right dimension.
            validN += 1;
        }

        //changing the sizes of the copies of A and B.
        copyA = allocateMatrix(validN,validN); 
        copyB = allocateMatrix(validN,validN);

        //filling the copies while putting zeros in the added lines and columns.
        for(int i = 0 ; i < validN ; i++){

          for(int j = 0 ; j < validN; j++){
            if((i < n) && (j < n)){  //the lines and columns of A and B
              copyA[(i*validN)+j] = A[(i*n) + j];
              copyB[(i*validN)+j] = B[(i*n) + j];

            }else{ //the new added lines and columns.
              copyA[(i*validN)+j] = 0;
              copyB[(i*validN)+j] = 0;
            }
          }
        }
        
    
        n = validN;
       
        

      }


      double* tmp,*tmp2;  //temporary variables used in strassen formulas.
      tmp = allocateMatrix(n/2,n/2);
      tmp2 = allocateMatrix(n/2,n/2);

      double* resBlock0,*resBlock1,*resBlock2,*resBlock3; //Elements of the result block.
      resBlock0 = allocateMatrix(n/2,n/2);
      resBlock1 = allocateMatrix(n/2,n/2);
      resBlock2 = allocateMatrix(n/2,n/2);
      resBlock3 = allocateMatrix(n/2,n/2);
      
      
      //Allocation for variables to stock the strassen formulas values.
      double* m1,*m2,*m3,*m4,*m5,*m6,*m7;
      m1 = allocateMatrix(n/2,n/2);
      m2 = allocateMatrix(n/2,n/2);
      m3 = allocateMatrix(n/2,n/2);
      m4 = allocateMatrix(n/2,n/2);
      m5 = allocateMatrix(n/2,n/2);
      m6 = allocateMatrix(n/2,n/2);
      m7 = allocateMatrix(n/2,n/2);

      
      //Variables to stock the smaller matrices.
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
       if(i < n/2){ //filling A0,A1,B0,B1.
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             A0[(i * (n/2)) + j] = copyA[(i * n) + j ] ;
             B0[(i * (n/2)) + j] = copyB[(i * n) + j ] ;
           }else{
             A1[(i * (n/2)) + (j % (n/2))] = copyA[(i * n) + j ] ;
             B1[(i * (n/2)) + (j % (n/2))] = copyB[(i * n) + j ] ;
           }
         }
       }else{ //filling A2,A3,B2,B3.
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             A2[((i%2) * (n/2)) + j] = copyA[(i * n) + j ] ;
             B2[((i%2) * (n/2)) + j] = copyB[(i * n) + j ] ;
           }else{
             A3[((i%2) * (n/2)) + (j % (n/2))] = copyA[(i * n) + j ] ;
             B3[((i%2) * (n/2)) + (j % (n/2))] = copyB[(i * n) + j ] ;
           }
         }
       }

        

      

      }

     
      //Calculating strassen formulas with the new matrices.
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


      //Calculating the result blocks.
      //resblock0
      matrixAdd(tmp,m1,m4,n/2,n/2);
      matrixSub(tmp2,m7,m5,n/2,n/2);
      matrixAdd(resBlock0,tmp,tmp2,n/2,n/2);
      //resblock1
      matrixAdd(resBlock1,m3,m5,n/2,n/2);
      //resblock2
      matrixAdd(resBlock2,m2,m4,n/2,n/2);
      //resblock3
      matrixAdd(tmp,m3,m6,n/2,n/2);
      matrixSub(tmp2,m1,m2,n/2,n/2);
      matrixAdd(resBlock3,tmp,tmp2,n/2,n/2);

      //filling the result matrix S from the blocks.
      for(int i = 0 ;i <= n-1; i++){
       if(i < n/2){ //filling resblock0,resblock1.
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             S[(i * n) + j ] =  resBlock0[(i * (n/2)) + j] ;
            
           }else{
             S[(i * n) + j ] = resBlock1[(i * (n/2)) + (j % (n/2))] ;
             
           }
         }
       }else{//filling resblock2,resblock3.
         for(int j = 0 ; j <= n-1 ; j++){
           if(j< n / 2){
             S[(i * n) + j ] = resBlock2[((i%2) * (n/2)) + j]  ;
             
           }else{
             S[(i * n) + j ] = resBlock3[((i%2) * (n/2)) + (j % (n/2))];
             
           }
         }
       }
      

      }

      //freeing the space dynamically allocated.
      free(copyA);
      free(copyB);

      free(tmp);
      free(tmp2);

      free(resBlock0);
      free(resBlock1);
      free(resBlock2);
      free(resBlock3);

      free(m1);
      free(m2);
      free(m3);
      free(m4);
      free(m5);
      free(m6);
      free(m7);

      free(A0);
      free(A1);
      free(A2);
      free(A3);
      free(B0);
      free(B1);
      free(B2);
      free(B3);

      


    }

    
}

/* 
    Solves a system of linear equations Ax=b for a double-precision matrix A (size n x n).
    Uses iterative ascension algorithm. 
    After the procedure, x contains the solution of Ax=b.
    We assume that x has been allocated outside the function.
*/
void        SolveTriangularSystemUP   (double *x, double *A, double *b, uint64_t n){

    for(int j = n-1 ; j >= 0 ; j-- ){ //index indicating the number of the column starting from the last column.
      
      uint64_t indexPivot = (j * n) + j;  //the pivot is in the diagonal of A so the index of the line is the same as the column.
      double pivot = A[indexPivot];  //the pivot of the line

      b[j] = b[j] / pivot;  //index of the column of the pivot can also be used for the line in the vector.
      x[j] = b[j]; //stocking the solution in x.
      A[indexPivot] = pivot / pivot; //to make the pivot equal to 1.

      for(int i = j-1 ; i >= 0; i--){ //index of the line starting from the line above the pivot.
        b[i] = b[i] - (A[i*n+j]*b[j]); //modifiyng the vector first because we need A[i*n+j].
        A[(i*n)+j] = A[(i*n)+j] - (A[(i*n)+j]*pivot);

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
  int toPivot = n +1; //the number we add to the index of the pivot to reach the next one.

  for(int i = 0 ; i <= n*n ; i = i + toPivot){  //running through the pivots.

    if(A[i] == 0){  //if a pivot is null, we need to exchange the line with another one that is not null in the same column.

      int j = i; //variable to stock the index of the correct line.

      while(A[j] == 0 && j < n*n){ //searching for a non-null line.
        j = j + n;
      }

      if(j > n*n){ //if no non-null lines have been found,
        return false; //there is no solution for the matrix.
      }

      for(int k = j%n ; k < n ; k++){  //swapping the lines
        double tmp = A[i + k] ;
        A[i+k] = A[j+k];
        A[j+k] = tmp;


      }

      int LinePivot = (int)(i / n); //the line of the supposed pivot.
      int indVect = (int)(j / n);  //index(in the vector) of the correct line of the pivot.

      //swapping the lines in the vector
      double tmp = b[LinePivot] ;
      b[LinePivot] = b[indVect];
      b[indVect] = tmp;

    }
  

    //Gauss Elimination
    for(int p = i+n; p <= n*n; p = p + n){ //running through the elements below the pivot. 

      double valUnderPivot = A[p]; //the value of the element.
      int lastInLine = n - (p%n); //variable to stop the (for) loop when we reach the end of a line.
      for(int y = 0; y < lastInLine ; y++){ //changing the values of the elements of the lines below the pivot.
        A[p+y] = A[p+y] - ((valUnderPivot/A[i])*A[i+y]);
      }

      uint64_t LinePivot = (int)(i / n); //the line of the supposed pivot.
      uint64_t indVect = (int)(p / n);  //index(in the vector) of the correct line of the pivot.

      b[indVect] = b[indVect] - (valUnderPivot/A[i])*b[LinePivot];  //applying the same operation to the vector

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
    //step 1: Triangularize
    cout<< "triangularizing the matrix..." << endl;
    bool triangular = Triangularize(A,b,n);
    cout<< "triangularized matrix :" << endl;
    
    if(triangular){
      cout<< "triangularized with success. triangularized matrix :" << endl;
      writeMatrix(stdout, A, n, n);
      //step 2: Solve the system
      SolveTriangularSystemUP(x,A,b,n);
      return true;

    }else{ // if there is no solution for the matrix.
      return false;
    }

    
}

bool        decompLU           (double *A, uint64_t n){
  int toPivot = n +1; //the number we add to the index of the pivot to reach the next one.

  for(int i = 0 ; i <= n*n ; i = i + toPivot){  //running through the pivots.

    if(A[i] == 0){  //if a pivot is null, we need to exchange the line with another one that is not null in the same column.

      int j = i; //variable to stock the index of the correct line.

      while(A[j] == 0 && j < n*n){ //searching for a non-null line.
        j = j + n;
      }

      if(j > n*n){ //if no non-null lines have been found,
        return false; //there is no solution for the matrix.
      }

      for(int k = j%n ; k < n ; k++){  //swapping the lines
        double tmp = A[i + k] ;
        A[i+k] = A[j+k];
        A[j+k] = tmp;


      }


    }
  

    //Gauss Elimination
    for(int p = i+n; p <= n*n; p = p + n){ //running through the elements below the pivot. 

      double valUnderPivot = A[p]; //the value of the element.
      int lastInLine = n - (p%n); //variable to stop the (for) loop when we reach the end of a line.
      A[p] = (valUnderPivot/A[i]);
      for(int y = 1; y < lastInLine ; y++){ //changing the values of the elements of the lines below the pivot.
        A[p+y] = A[p+y] - ((valUnderPivot/A[i])*A[i+y]);
      }
      
    }

  }
    return true;
}

double det(double *A, uint64_t n){
  //decompLU(A,n);
  writeMatrix(stdout,A,n,n);

  int toPivot = n +1;
  double det = 1;
  for(int i = 0 ; i < n*n ; i+= toPivot){

    det = det * A[i]; 
  }

  return det;
}

