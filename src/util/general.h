#ifndef UTIL_GENERAL_H_
#define UTIL_GENERAL_H_

#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>

#define MKL_Complex16 Complex
#include "tensor/complex.h"
#include "mkl.h"
#include "sys/time.h"
#include"omp.h"

#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

#define PI 3.14159265358979323846

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::istream;
using std::ostream;
using std::string;

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
// flag = 0 descending   flag = 1 ascending
// return array and index
template <typename T>
void QuickSort(T* array, int* index, int left, int right, int flag=0);

template <typename T>
void QuickSort(T* array, int left, int right, int flag=0);

void ReorderRelevantArray(int num_quantum, int num_table, int** quantum_table, 
        int* index);

//void RealTriMatrixDiag(double* );
//

void RealSymMatrixDiag(double* eigenvector, double* eigenvalue, int vector_dim);

void ComplexSymMatrixDiag(Complex* eigenvector, double* eigenvalue, int vector_dim);
// complex matrix diagonal
//
inline void CompMatrixDiagon(Complex* eigenVector,double* eigenValue, int vectorDimen)
{
    char Jobz,UpLo;
    Complex *ComWorkArea;
    double *DouWorkArea;
    int *IntWorkArea,ComWorkSize,DouWorkSize,IntWorkSize,Info;
    
    Jobz='V';UpLo='L';
    
    ComWorkSize=vectorDimen*vectorDimen+2*vectorDimen;
    DouWorkSize=2*vectorDimen*vectorDimen+5*vectorDimen+1;
    IntWorkSize=5*vectorDimen+3;
    
    ComWorkArea=new Complex[ComWorkSize];
    DouWorkArea=new double[DouWorkSize];
    IntWorkArea=new int[IntWorkSize];
    
    zheevd(&Jobz,&UpLo,&vectorDimen,eigenVector,&vectorDimen,eigenValue,ComWorkArea,&ComWorkSize,DouWorkArea,&DouWorkSize,IntWorkArea,&IntWorkSize,&Info);
    
    delete[] ComWorkArea;
    delete[] DouWorkArea;
    delete[] IntWorkArea;
}

inline double get_wall_time()
{
    struct timeval time;
    if(gettimeofday(&time,NULL)){
        return 0;
    }

    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

inline void error(const string& s)
{
    cout << endl << s << endl;
    cout.flush();
    exit(-1);
}

inline double GetWallTime()
{
    struct timeval time;
    if(gettimeofday(&time, NULL)) return 0;
    
    return (double)time.tv_sec + (double)time.tv_usec*0.000001;
}

inline bool CompareQuantumTable(int num_quantum, int* quantum_table1, int* quantum_table2)
{
    bool result = true;
    for(int i=0;i<num_quantum;++i)
    {
        if(quantum_table1[i] != quantum_table2[i])
        {
            result = false;
            break;
        }
    }
    return result;
}

inline void isParallelElement(double* element, int position, int& info, double& prefactor,
                              bool& zero)
{
    double factor;
    if(fabs(element[1]) <= 1E-10)
    {
        if(fabs(element[0]) > 1E-10)
        {
            info = -1;
            prefactor = 1.0;
        }
    }
    else
    {
        factor = element[0]/element[1];
        if(zero == true)
        {
            info = position;
            prefactor = factor;
            zero = false;
        }
        else if(zero == false)
        {
            if(fabs(factor-prefactor) > 1E-8)
            {
                info = -1;
                prefactor = 1.0;
            }
        }
    }
}

inline void isParallelElement(Complex* element, int position, int& info, Complex& prefactor,
                              bool& zero)
{
    Complex factor;
    if(Norm(element[1]) <= 1E-10)
    {
        if(Norm(element[0]) > 1E-10)
        {
            info = -1;
            prefactor = 1.0;
        }
    }
    else
    {
        factor = element[0]/element[1];
        if(zero == true)
        {
            info = position;
            prefactor = factor;
            zero = false;
        }
        else if(zero == false)
        {
            if(Norm(factor-prefactor) > 1E-8)
            {
                info = -1;
                prefactor = 1.0;
            }
        }
    }
}

inline void MkdirCacheFolder(bool disk_cache, char* cache_name)
{
    char label[2048];

    if(disk_cache == true)
    {
        if(cache_name[0] != '.')
        {
            sprintf(label, "mkdir %s", cache_name);
            system(label);
        }

        sprintf(label, "mkdir %s/SpaceDiskCacheFile", cache_name);
        system(label);
        sprintf(label, "mkdir %s/NetworkDiskCacheFile", cache_name);
        system(label);
    }
}


#endif // UTIL_GENERAL_H_

