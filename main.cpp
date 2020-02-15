#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "matrix.hpp"
#include <cstdint>
using namespace std;


int main(int argc, char** argv){
    /*
    //EX1
    double* A = allocateMatrix(2,2);
    A[0] = 0;
    A[1] = 10;
    A[2] = 0;
    A[3] = 0;

    double* B = allocateMatrix(2,2);
    B[0] = 0;
    B[1] = 10;
    B[2] = 0;
    B[3] = 0;
    */

    //EX2
    double* A = allocateMatrix(4,4);
    A[0] = 1;
    A[1] = 2;
    A[2] = 3;
    A[3] = 4;
    A[4] = 0;
    A[5] = 2;
    A[6] = 0;
    A[7] = 3;
    A[8] = 5;
    A[9] = 1;
    A[10] = 0;
    A[11] = 0;
    A[12] = 0;
    A[13] = 0;
    A[14] = 1;
    A[15] = 3;


    double* B = allocateMatrix(4,4);
    B[0] = 0;
    B[1] = 0;
    B[2] = 1;
    B[3] = 0;
    B[4] = 1;
    B[5] = 2;
    B[6] = 0;
    B[7] = 1;
    B[8] = 0;
    B[9] = 0;
    B[10] = 0;
    B[11] = 2;
    B[12] = 1;
    B[13] = 0;
    B[14] = 3;
    B[15] = 4;
    


    double* S = allocateMatrix(4,4);
    matrixMultiplyNaive(S,A,B,4,4,4);

    writeMatrix(stdout,S,4,4);

    

    return 0;
}