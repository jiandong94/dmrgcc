#ifndef DMRGCC_TENSOR_COMPLEX_MATRIX_BLOCK_H_
#define DMRGCC_TENSOR_COMPLEX_MATRIX_BLOCK_H_

#include "tensor/complex_matrix.h"

class ComplexMatrixBlock
{
    friend class ComplexTensorContraction;
    protected:

    int num_block_;
    
    int* left_index_;
    int* right_index_;

    ComplexMatrix** matrix_block_;


    public:

    // constructor
    //
    ComplexMatrixBlock();

    // constructor
    ComplexMatrixBlock(int num_block);

    // constructor
    //
    ComplexMatrixBlock(int num_block, int* left_index, int* right_index);

    // constructor
    //
    ComplexMatrixBlock(ComplexMatrixBlock* tmp_matrix_block);

    // destructor
    //
    ~ComplexMatrixBlock();

    // 
    //
    int get_num_block();

    //
    //
    int* get_left_index();

    //
    //
    int* get_right_index();

    //
    //
    ComplexMatrix* get_matrix_block(int position);

    //
    //
    void set_matrix_block(int position, ComplexMatrix* tmp_matrix);

    // compute the total dimension of matrices
    //
    int ComputeMatrixBlockDim();

    // compute the total dimension of matrices before position.
    // position = 0, return 0
    // position = 1, return dimension of the first matrix
    // position = 2, return dimension of the first two matrices
    int ComputePartMatrixBlockDim(int position);

    //
    //
    void NormalizeMatrixBlock();

    //
    // 
    void VectorizeMatrixBlock(bool direction, Complex* state);

    //
    //
    void PrintMatrixBlock();

    //
    //
    void WriteMatrixBlock(const char* matrix_block_name);

    //
    //
    void WriteMatrixBlock(ofstream &matrix_block_file);

    //
    //
    void ReadMatrixBlock(const char* matrix_block_name);

    //
    //
    void ReadMatrixBlock(ifstream &matrix_block_flie);

    //
    //
    void ResetMatrixBlock();

    //
    //
    void AddToMatrixBlock(int position, Complex factor, ComplexMatrix* tmp_matrix);

    //
    //
    void MultiplyToScalar(Complex scalar);

    // return position
    //
    int FindMatrixBlock(int left, int right);

    // find number and position of matrix with the same leigh target index
    //
    void FindMatrixBlock(int &num_same_index, int* &position_same_index, 
                         int leigh, int target_index);

};

#endif // DMRGCC_TENSOR_REAL_MATRIX_BLOCK_H_
