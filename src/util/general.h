#ifndef UTIL_GENERAL_H_
#define UTIL_GENERAL_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>

#include "mkl.h"
#include "sys/time.h"

#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

// Computes a matrix-matrix product with general matrices.
// C := alpha*op(A)*op(B) + beta*C.
// Layout = CblasRowMajor
// transa = CblasNoTrans
// transb = CblasNoTrans
// m = Specifies the number of rows of the matrix a
// n = Specifies the number of columns of the matrix b
// k = Specifies the number of columns of the matrix a
// alpha = Specifies the scalar alpha
// a = Array, size lda*m
// lda = lda must be at least max(1, k)
// b = Array, size ldb by k
// ldb = ldb must be at least max(1, n)
// c = Specifies the scalar beta. When beta is equal to zero, then c need not be set on input.
// ldc = ldc must be at least max(1, n)
//extern "C" void cblas_dgemm (CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb, 
//                            MKL_INT m, MKL_INT n, MKL_INT k, double alpha, double *a, 
//                             MKL_INT lda, double *b, MKL_INT ldb, double beta, double *c, MKL_INT ldc);

// quick sort
// return array and index
void QuickSort(double* array, int* index, int left, int right);


inline double GetWallTime()
{
    struct timeval time;
    if(gettimeofday(&time, NULL)) return 0;
    
    return (double)time.tv_sec + (double)time.tv_usec*0.000001;
}


#endif // UTIL_GENERAL_H_

