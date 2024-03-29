#ifndef DMRGCC_TENSOR_REAL_MATRIX_H_
#define DMRGCC_TENSOR_REAL_MATRIX_H_

//#include "util/general.h"
#include "util/general.ih"

//  A class of RealMatrix.
//  This class includes usual methods to operare real matrix.
//  The elements are stored in an array and row-major.
//  
class RealMatrix
{
    
    protected:

    int row_;
    int column_;
    int total_element_num_;

    double* matrix_element_;

    public:

    // constructor
    //
    RealMatrix();

    // constructor
    //
    RealMatrix(int row, int column);
    
    // copy constructor
    //
    RealMatrix(RealMatrix* tmp_matrix);

    // destructor
    //
    ~RealMatrix();

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
    void set_matrix_element(int row, int column, double element);

    // get a matrix element
    // 
    double get_matrix_element(int row, int column);
    
    // get all the matrix elements
    // 
    double* get_matrix_element();


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
    void AddToMatrixElement(int row, int column, double element);

    // add the matrix to matrix
    //
    void AddToMatrix(RealMatrix* tmp_matrix);
   
    // matrix + factor*tmp_matrix
    //
    void AddToMatrix(double factor, RealMatrix* tmp_matrix);

    // matrix * scalar
    //
    void MultiplyToScalar(double scalar);

    // matrix * tmp_matrix
    //
    RealMatrix* MultiplyToMatrix(RealMatrix* tmp_matrix);

    // matrix (kronecker product) tmp_matrix
    //
    RealMatrix* MatrixKronProduct(RealMatrix* tmp_matrix);

    // matrix_element * matrix_element
    //
    void MatrixElementProduct(RealMatrix* tmp_matrix);

    // change matrix
    // if leigh = 0, change matrix dimension to (dimen, RightDimen)
    //     if dimen >= LeftDimen, put the matrix elements in the new matrix,
    //     else change the new matrix elements to zero.
    void ChangeMatrix(int leight, int truncate_dim);

    // matrix'
    //
    RealMatrix* TransposeMatrix();

    // reshape matrix to (row, column)
    //
    RealMatrix* ReshapeMatrix(int row, int column);

    // replace matrix(copy matrix)
    //
    void ReplaceMatrix(RealMatrix* tmp_matrix);
    
    // deparallelisation algorithm
    // M = M'*T (T is transfer matrix)
    // M' has at most as many columns as M and no two columns which are parallel to each other
    // eg:
    // [1 3 1;     [1 3;   [1 0 1;
    //          =        * 
    //  2 4 2]      2 4]    0 0 0]
    void ParallelMatrix(int leigh, int &num_unparallel, int* &position_unparallel, 
            RealMatrix* &transfer_tensor);

    // expan matrix
    // if flag = 0, expan the matrix to (row_+row, column_). 
    // i.e. piece the up-matrix and the down-matrix
    //
    void ExpanMatrix(int flag, RealMatrix* tmp_matrix);

    // svd
    void SVDMatrix(RealMatrix* &left_matrix, RealMatrix* &right_matrix, 
                   double* &singular_value, int &singular_dim);

};

#endif // DMRGCC_TENSOR_REAL_MATRIX_H_
