#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "matrix.hpp"
#include <cstdint>
using namespace std;


int main(int argc, char** argv){
    /*
    cout << "Test: matrixMultiplyNaive | expected result: S = { 9 , 2 , 33 , 4}" << endl;
    double* A = allocateMatrix(2,2);
    A[0] = 1;
    A[1] = 2;
    A[2] = 5;
    A[3] = 4;

    double* B = allocateMatrix(2,2);
    B[0] = 5;
    B[1] = 0;
    B[2] = 2;
    B[3] = 1;

    double* S = allocateMatrix(2,2);
    matrixMultiplyNaive(S,A,B,2,2,2);
    writeMatrix(stdout, S, 2, 2);

    */
   
   /*
    cout << "Test: matrixMultiplyNaive | expected result: S = { 7 , 6 , 9 , 11 }" << endl;
    double* A = allocateMatrix(2,3);
    A[0] = 1;
    A[1] = 0;
    A[2] = 2;
    A[3] = 3;
    A[4] = 1;
    A[5] = 2;

    double* B = allocateMatrix(3,2);
    B[0] = 1;
    B[1] = 2;
    B[2] = 0;
    B[3] = 1;
    B[4] = 3;
    B[5] = 2;

    double* S = allocateMatrix(2,2);
    matrixMultiplyNaive(S,A,B,2,3,2);
    writeMatrix(stdout, S, 2, 2);
    
 */

    /*

    cout << "Test: matrixMultiplyStrassen with 2*2 matrices |  expected result: S = { 9 , 2 , 33 , 4}" << endl;
    double* A = allocateMatrix(2,2);
    A[0] = 1;
    A[1] = 2;
    A[2] = 5;
    A[3] = 4;

    double* B = allocateMatrix(2,2);
    B[0] = 5;
    B[1] = 0;
    B[2] = 2;
    B[3] = 1;

    double* S = allocateMatrix(2,2);
    matrixMultiplyStrassen(S,A,B,2);
    writeMatrix(stdout, S, 2, 2);

    */

   /*

    cout << "Test: matrixMultiplyStrassen with 3*3 matrices |  expected result: S = {4,0,5,7,0,5,8,1,14}" << endl;
    double* A = allocateMatrix(3,3);
    A[0] = 1;
    A[1] = 0;
    A[2] = 2;
    A[3] = 3;
    A[4] = 0;
    A[5] = 1;
    A[6] = 2;
    A[7] = 1;
    A[8] = 4;

    double* B = allocateMatrix(3,3);
    B[0] = 2;
    B[1] = 0;
    B[2] = 1;
    B[3] = 0;
    B[4] = 1;
    B[5] = 4;
    B[6] = 1;
    B[7] = 0;
    B[8] = 2;

    double* S = allocateMatrix(4,4);
    matrixMultiplyStrassen(S,A,B,3);
    writeMatrix(stdout, S, 4, 4);
    
    */

    cout << "Test: SolveTriangularSystemUP | expected result: {-97/200 , -11/20 , 3/4 , 1/2}  "<< endl;

        double* A = allocateMatrix(4,4);
        A[0] = 10;
        A[1] = 3;
        A[2] = 8;
        A[3] = 7;
        A[4] = 0;
        A[5] = 5;
        A[6] = 3;
        A[7] = 5;
        A[8] = 0;
        A[9] = 0;
        A[10] = 2;
        A[11] = 3;
        A[12] = 0;
        A[13] = 0;
        A[14] = 0;
        A[15] = 6;

        double* x = allocateMatrix(4,1);
        setMatrixZero(x,4,1);

        double* b = allocateMatrix(4,1);
        b[0] = 3;
        b[1] = 2;
        b[2] = 3;
        b[3] = 3;

        SolveTriangularSystemUP(x,A,b,4);
        writeMatrix(stdout, x, 4, 1);







    /*
    //EX1
    double* A = allocateMatrix(2,2);
    A[0] = 2;
    A[1] = 3;
    A[2] = 1;
    A[3] = 4;

    double* B = allocateMatrix(2,2);
    B[0] = 0;
    B[1] = 2;
    B[2] = 3;
    B[3] = 0;
    */
    /*
    //EX2
    double* A = allocateMatrix(4,4);
    A[0] = 1.0;
    A[1] = 2.0;
    A[2] = 3.0;
    A[3] = 4.0;
    A[4] = 0.0;
    A[5] = 2.0;
    A[6] = 0.0;
    A[7] = 3.0;
    A[8] = 5.0;
    A[9] = 1.0;
    A[10] = 0.0;
    A[11] = 0.0;
    A[12] = 0.0;
    A[13] = 0.0;
    A[14] = 1.0;
    A[15] = 3.0;

    /*
    double* B = allocateMatrix(4,4);
    B[0] = 0.0;
    B[1] = 0.0;
    B[2] = 1.0;
    B[3] = 0.0;
    B[4] = 1.0;
    B[5] = 2.0;
    B[6] = 0.0;
    B[7] = 1.0;
    B[8] = 0.0;
    B[9] = 0.0;
    B[10] = 0.0;
    B[11] = 2.0;
    B[12] = 1.0;
    B[13] = 0.0;
    B[14] = 3.0;
    B[15] = 4.0;
    */
    /*
    //EX3
    double* A = allocateMatrix(3,2);
    A[0] = 2;
    A[1] = 4;
    A[2] = -2;
    A[3] = 4;
    A[4] = -2;
    A[5] = 6;
    
    
   
    

    double* B = allocateMatrix(2,3);
    B[0] = 6;
    B[1] = -4;
    B[2] = 2;
    B[3] = 2;
    B[4] = 5;
    B[5] = 4;
    */

    /*
    double* S = allocateMatrix(3,3);
    matrixMultiplyNaive(S,A,B,3,2,3);

    writeMatrix(stdout,S,3,3);
    */
   /*
    cout << "Test: matrixMultiplyNaive" << endl;
    cout << "expected result" << endl; 
    double res [] = {-2,-1,10,17,7,-22,-17,-7,34};
    writeMatrix(stdout, res, 3, 3);
    double A [] = {3,1,-2,-1,0,4,5,2,-6};
    double B [] = {3,1,-2,-1,0,4,5,2,-6};
    double S [] = {0,0,0,0,0,0,0,0,0};
    matrixMultiplyNaive(S,A,B,3,3,3);
    cout << "obtained result:" << endl; 
    writeMatrix(stdout, S, 3, 3);
    */
    /*
    cout << "Test: Triangularize | expected result: A = [1,1,0,1], b = [1,1]" << endl;
    double A [] = {4,1,1,1};
    double b [] = {1,1};
    Triangularize(A,b,2);
    cout << "The matrix A : " << endl;
    writeMatrix(stdout, A, 2,2);
    cout << "The vector b : " << endl;
    writeMatrix(stdout, b, 2, 1);
    */
    /*
    cout << "Test: Triangularize | expected result: A[1,0,1,0,3,1,0,0,0.333] b[1,1,0.333]" << endl;
    double A [] = {0,2,1,0,3,1,1,0,1};
    double b [] = {1,1,1};
    Triangularize(A,b,3);
    cout << "The matrix A : " << endl;
    writeMatrix(stdout, A, 3,3);
    cout << "The vector b : " << endl;
    writeMatrix(stdout, b, 3, 1);
    */

    /*    
    cout << "Test: SolveTriangularSystemUP | expected result: 1 , 2 , 3" << endl;
        double A [] = {0,5,4,8,2,0,8,0,3};
        double x [] = {0,0,0};
        double b [] = {1,4,9};
        SolveSystemGauss(x,A,b,3);
        writeMatrix(stdout, x, 1, 3);
    */
   /*
   cout << "Test: solvouhjokeStrassen" << endl;
        double* A = allocateMatrix(3,3);
        A[0] = 1;
        A[1] = 0;
        A[2] = 1;
        A[3] = 5;
        A[4] = 0;
        A[5] = 4;
        A[6] = 2;
        A[7] = 1;
        A[8] = 0;
        double* B = allocateMatrix(3,3);
        B[0] = 0;
        B[1] = 1;
        B[2] = 2;
        B[3] = 1;
        B[4] = 2;
        B[5] = 4;
        B[6] = 0;
        B[7] = 2;
        B[8] = 1;
        
        
        double* S = allocateMatrix(4,4);
        matrixMultiplyStrassen(S,A,B,3);
        writeMatrix(stdout, S, 4, 4);
       

    */
    return 0;
}