#include "tensor/real_matrix_block.h"

int main()
{
    cout << "=================================" << endl;
    cout << "       Test RealMatrixBlock      " << endl;
    cout << "=================================" << endl;

    //srand((unsigned)time(NULL));


    cout << "1. Constructor" << endl;
    // constructor
    RealMatrixBlock* matrix_block;
    matrix_block = new RealMatrixBlock();
    matrix_block->PrintMatrixBlock();
    delete matrix_block;

    // constructor
    matrix_block = new RealMatrixBlock(4);
    matrix_block->PrintMatrixBlock();
   
    RealMatrix* tmp_matrix = new RealMatrix(3,3);
    for(int i=0;i<matrix_block->get_num_block();++i)
    {
        tmp_matrix->RandomMatrix();
        matrix_block->set_matrix_block(i, tmp_matrix);
    }
    matrix_block->PrintMatrixBlock();
    delete tmp_matrix;


    cout << endl;
    cout << "2. Get" << endl;
    int num_block;
    int *left_index, *right_index;
    RealMatrix* matrix_0;
    matrix_block->PrintMatrixBlock();
    num_block = matrix_block->get_num_block();
    left_index = matrix_block->get_left_index();
    right_index = matrix_block->get_right_index();
    matrix_0 = matrix_block->get_matrix_block(0);
    cout << "number of block: " << num_block << endl;
    cout << "left index: " ;
    for(int i=0;i<num_block;++i) cout << left_index[i] << " ";
    cout << endl;
    cout << "right index: " ;
    for(int i=0;i<num_block;++i) cout << right_index[i] << " ";
    cout << endl;
    cout << "the first matrix: " << endl;
    matrix_0->PrintMatrix();


    cout << endl;
    cout << "3. ComputeMatrixBlockDim ComputePartMatrixBlockDim " << endl;
    cout << "   NormalizeMatrixBlock  VectorizeMatrixBlock" << endl;
    matrix_block->PrintMatrixBlock();
    cout << "compute matrix block dim: " << matrix_block->ComputeMatrixBlockDim() << endl;
    cout << "compute part matrix block dim: first matrix  " << matrix_block->ComputePartMatrixBlockDim(1) << endl;
    
    cout << "NormalizeMatrixBlock" << endl;
    matrix_block->NormalizeMatrixBlock();
    matrix_block->PrintMatrixBlock();
    double sumsquare=0;
    for(int i=0;i<matrix_block->get_num_block();++i)
    {
        tmp_matrix = matrix_block->get_matrix_block(i);
        sumsquare += tmp_matrix->SumSquareMatrix();
    }
    cout << "sum square: " << sumsquare << endl;

    cout << "VectorizeMatrixBlock" << endl;
    double* state = new double[matrix_block->ComputeMatrixBlockDim()];
    matrix_block->VectorizeMatrixBlock(true, state);
    cout << "matrix_block->state: " << endl;
    for(int i=0;i<matrix_block->ComputeMatrixBlockDim();++i) cout << state[i] << " " << endl;

    for(int i=0;i<matrix_block->ComputeMatrixBlockDim();++i) state[i] = i+1;
    cout << "state->matrix_block: " << endl;
    matrix_block->VectorizeMatrixBlock(false, state);
    matrix_block->PrintMatrixBlock();

    delete[] state;


    cout << endl;
    cout << "4. WriteMatrixBlock ReadMatrixBlock" << endl;
    char const *char_matrix_block = "matrix_block.dat";
    matrix_block->PrintMatrixBlock();
    cout << "Write matrix block ..." << endl;
    matrix_block->WriteMatrixBlock(char_matrix_block);
    cout << "Read matrix block ..." << endl;
    RealMatrixBlock* matrix_block_read = new RealMatrixBlock();
    matrix_block_read->ReadMatrixBlock(char_matrix_block);
    matrix_block_read->PrintMatrixBlock();
    delete matrix_block_read;



    cout << endl;
    cout << "5. FindMatrixBlock" << endl;
    cout << "position of index (1,1): " <<matrix_block->FindMatrixBlock(1,1) << endl;

    int num_same_index;
    int* position_same_index;
    matrix_block->FindMatrixBlock(num_same_index, position_same_index, 0, 2);
    cout << "number of same index for left index 2: " << num_same_index << endl;
    cout << "position of same index for left index 2: " << position_same_index[0] << endl;

    delete[] position_same_index;
    
    delete matrix_block;
    return 0;

}
