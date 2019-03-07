#ifndef DMRGCC_TENSOR_COMPLEX_MATRIX_H_
#define DMRGCC_TENSOR_COMPLEX_MATRIX_H_

#include "tensor/complex.h"
#include "util/general.ih"

//  A class of ComplexMatrix.
//  This class includes usual methods to operate complex matrix.
//  The elements are stored in an array and row-major.
//  
class ComplexMatrix
{
    
    protected:

    int row_;
    int column_;
    int total_element_num_;

    Complex* matrix_element_;

    public:

    // constructor
    //
    ComplexMatrix();

    // constructor
    //
    ComplexMatrix(int row, int column);
    
    // copy constructor
    //
    ComplexMatrix(ComplexMatrix* tmp_matrix);

    // destructor
    //
    ~ComplexMatrix();

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
    void set_matrix_element(int row, int column, Complex element);

    // get a matrix element
    // 
    Complex get_matrix_element(int row, int column);
    
    // get all the matrix elements
    // 
    Complex* get_matrix_element();


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

    // sum of square of all matrix elements
    //
    double SumSquareMatrix();
    
    // add a value to a matrix element
    //
    void AddToMatrixElement(int row, int column, Complex element);

    // add the matrix to matrix
    //
    void AddToMatrix(const ComplexMatrix* tmp_matrix);
   
    // matrix + factor*tmp_matrix
    //
    void AddToMatrix(Complex factor, const ComplexMatrix* tmp_matrix);

    // matrix * scalar
    //
    void MultiplyToScalar(Complex scalar);

    // matrix * tmp_matrix
    //
    ComplexMatrix* MultiplyToMatrix(ComplexMatrix* tmp_matrix);

    // matrix (kronecker product) tmp_matrix
    //
    ComplexMatrix* MatrixKronProduct(ComplexMatrix* tmp_matrix);

    // matrix_element * matrix_element
    //
    void MatrixElementProduct(ComplexMatrix* tmp_matrix);

    // change matrix
    // if leigh = 0, change matrix dimension to (dimen, RightDimen)
    //     if dimen >= LeftDimen, put the matrix elements in the new matrix,
    //     else change the new matrix elements to zero.
    void ChangeMatrix(int leight, int truncate_dim);

    // matrix'
    //
    ComplexMatrix* TransposeMatrix();

    //
    //
    ComplexMatrix* HermitianConjugateMatrix();

    // reshape matrix to (row, column)
    //
    ComplexMatrix* ReshapeMatrix(int row, int column);

    // replace matrix(copy matrix)
    //
    void ReplaceMatrix(ComplexMatrix* tmp_matrix);
    
    // deparallelisation algorithm
    // M = M'*T (T is transfer matrix)
    // M' has at most as many columns as M and no two columns which are parallel to each other
    // eg:
    // [1 3 1;     [1 3;   [1 0 1;
    //          =        * 
    //  2 4 2]      2 4]    0 0 0]
    void ParallelMatrix(int leigh, int &num_unparallel, int* &position_unparallel, 
            ComplexMatrix* &transfer_tensor);

    // expan matrix
    // if flag = 0, expan the matrix to (row_+row, column_). 
    // i.e. piece the up-matrix and the down-matrix
    //
    void ExpanMatrix(int flag, ComplexMatrix* tmp_matrix);

    // svd
    void SVDMatrix(ComplexMatrix* &left_matrix, ComplexMatrix* &right_matrix, 
                   double* &singular_value, int &singular_dim);

};

#endif // DMRGCC_TENSOR_REAL_MATRIX_H_
