#include "tensor/real_matrix_block.h"

int main()
{
    srand((unsigned)time(NULL));

    // constructor
    RealMatrixBlock* matrix_block1;
    matrix_block1 = new RealMatrixBlock();
    delete matrix_block1;

    // constructor
    matrix_block1 = new RealMatrixBlock(4);
    cout << "number of matrix block: " << matrix_block1->get_num_block() << endl;
    
    // Reset
    matrix_block1->ResetMatrixBlock();
    RealMatrix* tmp_matrix = new RealMatrix(3,3);
    for(int i=0;i<matrix_block1->get_num_block();++i)
    {
        tmp_matrix->RandomMatrix();
        matrix_block1->set_matrix_block(i, tmp_matrix);
    }

    matrix_block1->PrintMatrixBlock();
    delete tmp_matrix;


    // find matrix
    cout << "position of index (1,1): " <<matrix_block1->FindMatrixBlock(1,1) << endl;

    int num_same_index;
    int* position_same_index;
    matrix_block1->FindMatrixBlock(num_same_index, position_same_index, 0, 2);
    cout << "number of same index for left index 2: " << num_same_index << endl;
    cout << "position of same index for left index 2: " << position_same_index[0] << endl;

    delete[] position_same_index;
    
    // read write
    cout << "read and write matrix block" << endl;
    char* test_file = "matrix_block.dat";
    matrix_block1->PrintMatrixBlock();
    matrix_block1->WriteMatrixBlock(test_file);

    RealMatrixBlock* matrix_block2 = new RealMatrixBlock();
    matrix_block2->ReadMatrixBlock(test_file);
    matrix_block2->PrintMatrixBlock();
    
    
    delete matrix_block1;
    return 0;




}
