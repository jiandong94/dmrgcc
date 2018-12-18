#ifndef DMRGCC_TENSOR_REAL_MATRIX_BLOCK_H_
#define DMRGCC_TENSOR_REAL_MATRIX_BLOCK_H_

#include "tensor/real_matrix.h"

class RealMatrixBlock
{
    protected:

    int num_block_;
    
    int* left_index_;
    int* right_index_;

    RealMatrix** matrix_block_;


    public:

    // constructor
    //
    RealMatrixBlock();

    // constructor
    RealMatrixBlock(int num_block);

    // constructor
    //
    RealMatrixBlock(int num_block, int* left_index, int* right_index);

    // constructor
    //
    RealMatrixBlock(RealMatrixBlock* tmp_matrix_block);

    // destructor
    //
    ~RealMatrixBlock();

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
    RealMatrix* get_matrix_block(int position);

    //
    //
    void set_matrix_block(int position, RealMatrix* tmp_matrix);

    //
    //
    void PrintMatrixBlock();

    //
    //
    void WriteMatrixBlock(ofstream &matrix_block_file);

    //
    //
    void WriteMatrixBlock(char* matrix_block_name);

    //
    //
    void ReadMatrixBlock(ifstream &matrix_block_flie);

    //
    //
    void ReadMatrixBlock(char* matrix_block_name);

    //
    //
    void ResetMatrixBlock();

    // return position
    //
    int FindMatrixBlock(int left, int right);

    //
    //
    void FindMatrixBlock(int &num_same_index, int* &position_same_index, 
                         int leigh, int target_index);




};

#endif // DMRGCC_TENSOR_REAL_MATRIX_BLOCK_H_
