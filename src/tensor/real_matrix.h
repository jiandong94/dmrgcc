#ifndef DMRGCC_TENSOR_REAL_MATRIX_H_
#define DMRGCC_TENSOR_REAL_MATRIX_H_

#include "util/general.h"

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
    void WriteMatrix(char* matrix_name);

    // write
    //
    void WriteMatrix(ofstream &matrix_file);

    // read
    //
    void ReadMatrix(char* matrix_name);

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
    
    // add a value to a matrix element
    //
    void AddToMatrixElement(int row, int column, double element);

    // add the matrix to matrix
    //
    void AddToMatrix(RealMatrix* tmp_matrix);

    // multiply the matrix to scalar
    //
    void MultiplyToScalar(double scalar);

    // multiply the matrix to matrix
    // the real matrix product
    //
    RealMatrix* MultiplyToMatrix(RealMatrix* tmp_matrix);

    // element-wise product
    //
    void MatrixElementProduct(RealMatrix* tmp_matrix);

    //
    //
    RealMatrix* TransposeMatrix();

    // reshape matrix to (row, column)
    //
    RealMatrix* ReshapeMatrix(int row, int column);

    // expan matrix
    // if flag = 0, expan the matrix to (row_+row, column_). 
    // i.e. piece the up-matrix and the down-matrix
    //
    void CombineMatrix(int flag, RealMatrix* tmp_matrix);

    // svd
    void SVDMatrix(RealMatrix* &left_matrix, RealMatrix* &right_matrix, 
                   double* &sigular_value, int &sigular_dim);


};

#endif // DMRGCC_TENSOR_REAL_MATRIX_H_
