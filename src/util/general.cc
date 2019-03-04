#include "util/general.h"
#include "util/general.ih"

void ReorderRelevantArray(int num_quantum, int num_table, 
        int** quantum_table, int* index)
{
    int **tmp_quantum_table;
    tmp_quantum_table = new int* [num_table];
    for(int i=0;i<num_table;++i)
        tmp_quantum_table[i] = new int[num_quantum];
    // reorder
    int *position;
    position = new int[num_table];
    for(int i=0;i<num_table;++i)
    {
        position[index[i]] = i; // (index, position)
    }
    for(int i=0;i<num_table;++i)
        for(int j=0;j<num_quantum;++j)
            tmp_quantum_table[position[i]][j] = quantum_table[i][j];
    /*
    for(int i=0;i<num_table;++i)
        for(int j=0;j<num_quantum;++j)
            for(int k=0;k<num_table;++k)
            {
                if(index[k] == i)
                    tmp_quantum_table[k][j] = quantum_table[i][j];
            }
    */
    // quantum_table = tmp_quantum_table
    for(int i=0;i<num_table;++i)
        for(int j=0;j<num_quantum;++j)
            quantum_table[i][j] = tmp_quantum_table[i][j];

    for(int i=0;i<num_table;++i) delete[] tmp_quantum_table[i];
    delete[] tmp_quantum_table;
    delete[] position;
}

//void RealTriMatrixDiag()
//{
//    int matrix_layout = LAPACK_ROW_MAJOR; 
//    char jobz = 'V';
//}


void RealSymMatrixDiag(double* eigenvector, double* eigenvalue, int vector_dim)
{
    int matrix_layout = LAPACK_ROW_MAJOR;
    char jobz, uplo;
    lapack_int n, lda, info;

    jobz = 'V';
    uplo = 'U';
    n = vector_dim;
    lda = n;

    info = LAPACKE_dsyevd (matrix_layout, jobz, uplo, n, eigenvector,
            lda, eigenvalue);
    if(info != 0)
    {
        cout << "info = " << info << " in RealSymMatrixDiag" << endl;
    }
}

void ComplexSymMatrixDiag(Complex* eigenvector, double* eigenvalue, int vector_dim)
{
    int matrix_layout = LAPACK_ROW_MAJOR;
    char jobz, uplo;
    lapack_int n, lda, info;

    jobz = 'V';
    uplo = 'U';
    n = vector_dim;
    lda = n;

    info = LAPACKE_zheevd (matrix_layout, jobz, uplo, n, eigenvector,
            lda, eigenvalue);
    if(info != 0)
    {
        cout << "info = " << info << " in ComplexSymMatrixDiag" << endl;
    }
}
