#ifndef DMRGCC_TENSOR_MATRIX_H_
#define DMRGCC_TENSOR_MATRIX_H_

//#include "util/general.h"
#include "util/general.ih"
#include "tensor/complex.h"

//  A class of RealMatrix.
//  This class includes usual methods to operare real matrix.
//  The elements are stored in an array and row-major.
//  
template<class T>
class Matrix
{
    
    protected:

    int row_;
    int column_;
    int total_element_num_;

    T* matrix_element_;

    public:

    // constructor
    //
    Matrix();

    // constructor
    //
    Matrix(int row, int column);
    
    // copy constructor
    //
    Matrix(Matrix* tmp_matrix);

    // destructor
    //
    ~Matrix();

    // get row number
    //
    int get_row();

    // get column number
    //
    int get_column();

    // get total element number
    //
    int get_total_element_num();

    // set matrix element
    //
    void set_matrix_element(int row, int column, T element);

    // get a matrix element
    // 
    T get_matrix_element(int row, int column);
    
    // get all the matrix elements
    // 
    T* get_matrix_element();

    // print it
    //
    void PrintMatrix();

    // write
    //
    void WriteMatrix(const char* matrix_name);

    // write
    //
    void WriteMatrix(ofstream &matrix_file);

    // read
    //
    void ReadMatrix(const char* matrix_name);

    // read
    //
    void ReadMatrix(ifstream &matrix_file);


    // initial the whole matrix
    // set row and column to zero and matrix to ptrnull
    //
    void ResetMatrix();

    // reset all the matrix elements to zero
    //
    void ClearMatrix();

    // set matrix elements random
    //
    void RandomMatrix();

    //
    //
    double SumSquareMatrix();
    
    // add a value to a matrix element
    //
    void AddToMatrixElement(int row, int column, T element);

    // add the matrix to matrix
    //
    void AddToMatrix(Matrix* tmp_matrix);
    
    void AddToMatrix(T factor, Matrix* tmp_matrix);

    // multiply the matrix to scalar
    //
    void MultiplyToScalar(T scalar);

    // multiply the matrix to matrix
    // the real matrix product
    //
    Matrix* MultiplyToMatrix(Matrix* tmp_matrix);

    //
    //
    Matrix* MatrixKronProduct(Matrix* tmp_matrix);

    // element-wise product
    //
    void MatrixElementProduct(Matrix* tmp_matrix);

    // change matrix
    // if leigh = 0, change matrix dimension to (dimen, RightDImen)
    //     if dimen >= LeftDimen, put the matrix elements in the new matrix,
    //     else change the new matrix elements to zero.
    void ChangeMatrix(int leight, int truncate_dim);

    //
    //
    Matrix* TransposeMatrix();

    //
    //
    Matrix* HerimitMatrix();

    // reshape matrix to (row, column)
    //
    Matrix* ReshapeMatrix(int row, int column);

    // replace matrix(copy matrix)
    //
    void ReplaceMatrix(Matrix* tmp_matrix);
    
    // deparallelisation algorithm
    //
    void ParallelMatrix(int leigh, int &num_unparallel, int* &position_unparallel, 
            Matrix* &transfer_tensor);

    // expan matrix
    // if flag = 0, expan the matrix to (row_+row, column_). 
    // i.e. piece the up-matrix and the down-matrix
    //
    void ExpanMatrix(int flag, Matrix* tmp_matrix);

    // svd
    void SVDMatrix(Matrix* &left_matrix, Matrix* &right_matrix, 
                   double* &singular_value, int &singular_dim);


};

#endif // DMRGCC_TENSOR_REAL_MATRIX_H_
