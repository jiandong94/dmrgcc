#include "dmrg/real_tensor_lattice.h"

int main()
{
    cout << "=================================" << endl;
    cout << "       Test RealTensorLattice    " << endl;
    cout << "=================================" << endl;
    
    cout << "1. Constructor" << endl;
    RealTensorLattice* tensor_lattice = new RealTensorLattice(2);
    tensor_lattice->PrintTensorLattice();

    cout << endl;
    cout << "DefineTensorLattice" << endl;
    int num_left_block = 2;
    int num_right_block = 1;
    int *left_block = new int[num_left_block];
    int *right_block = new int[num_right_block];
    int *left_dim = new int[num_left_block];
    int *right_dim = new int[num_right_block];
    for(int i=0;i<num_left_block;++i)
    {
        left_block[i] = i;
        right_block[i] = i;
        left_dim[i] = 3;
        right_dim[i] = 1;
    }
    tensor_lattice->DefineTensorLattice(num_left_block, num_right_block, 
            left_block, right_block, left_dim, right_dim);
    tensor_lattice->PrintTensorLattice();
    delete[] left_block, right_block, left_dim, right_dim;

    cout << endl;
    cout << "DefineTensorLattice (ket tensor)" << endl;
    //int num_block = 4, row = 2, column = 3;
    //int left_index[] = {0, 1, 1, 2};
    //int right_index[] = {0, 1, 1, 2};
    //int physics_index[] = {0, 0, 0, 0};
    int num_block = 2, row = 3, column = 1;
    int left_index[] = {0, 1};
    int right_index[] = {0, 0};
    int physics_index[] = {0, 1};
    tensor_lattice->DefineTensorLattice(num_block, left_index, right_index, 
            physics_index, row, column);
    tensor_lattice->PrintTensorLattice();

    cout << endl;
    cout << "2. Get" << endl;
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();
    left_block = tensor_lattice->get_left_block();
    right_block = tensor_lattice->get_right_block();
    left_dim = tensor_lattice->get_left_dim();
    right_dim = tensor_lattice->get_right_dim();
    int **phy_index = tensor_lattice->get_physics_index();
    int physics_dim = tensor_lattice->get_physics_dim();
    cout << "physics_dim: " << physics_dim << endl;
    cout << "num_left_block: " << num_left_block << endl;
    cout << "num_right_block: " << num_right_block << endl;
    cout << "left_block: ";
    for(int i=0;i<num_left_block;++i) cout << left_block[i] << " ";
    cout << endl;
    cout << "right_block: ";
    for(int i=0;i<num_right_block;++i) cout << right_block[i] << " ";
    cout << endl;
    cout << "left_dim: ";
    for(int i=0;i<num_left_block;++i) cout << left_dim[i] << " ";
    cout << endl;
    cout << "right_block: ";
    for(int i=0;i<num_right_block;++i) cout << right_dim[i] << " ";
    cout << endl;
    cout << "physics_index:" << endl;
    for(int i=0;i<num_left_block;++i) 
    {
        for(int j=0;j<num_right_block;++j)
            cout << phy_index[i][j] << " ";
        cout << endl;
    }
       

    cout << endl;
    cout << "3. WriteTensorLattice ReadTensorLattice" << endl;
    char const *char_tensor_lattice = "tensor_lattice.dat";
    tensor_lattice->PrintTensorLattice();
    cout << "Write tensor lattice ..." << endl;
    tensor_lattice->WriteTensorLattice(char_tensor_lattice);
    cout << "Read tensor lattice ..." << endl;
    RealTensorLattice* tensor_lattice_read = new RealTensorLattice();
    tensor_lattice_read->ReadTensorLattice(char_tensor_lattice);
    tensor_lattice_read->PrintTensorLattice();
    delete tensor_lattice_read;
   
    cout << endl;
    cout << "4. ComputeLatticeDim ComputeKetTensorDim ComputePartKetTensorDim " << endl;
    cout << "   NormalizeTensorLattice  VectorizeTensorLattice" << endl;
    tensor_lattice->PrintTensorLattice();
    cout << "compute left bond dim: " << tensor_lattice->ComputeLatticeDim(0) << endl;
    cout << "compute right bond dim: " << tensor_lattice->ComputeLatticeDim(1) << endl;
    cout << "compute ket tensor dim: " << tensor_lattice->ComputeKetTensorDim() << endl;
    cout << "compute part ket tensor dim: first matrix  " << tensor_lattice->ComputePartKetTensorDim(1) << endl;

    cout << "NormalizeTensorLattice" << endl;
    tensor_lattice->NormalizeTensorLattice();
    tensor_lattice->PrintTensorLattice();
    //double sumsquare=0;
    //for(int i=0;i<tensor_lattice->get_ket_tensor()->matrix_block->get_num_block();++i)
    //{
    //   tmp_matrix = matrix_block->get_matrix_block(i);
    //    sumsquare += tmp_matrix->SumSquareMatrix();
    //}
    //cout << "sum square: " << sumsquare << endl;

    cout << "VectorizeTensorLattice" << endl;
    double* state = new double[tensor_lattice->ComputeKetTensorDim()];
    tensor_lattice->VectorizeTensorLattice(true, state);
    cout << "ket_tensor->state: " << endl;
    for(int i=0;i<tensor_lattice->ComputeKetTensorDim();++i) cout << state[i] << " " << endl;

    for(int i=0;i<tensor_lattice->ComputeKetTensorDim();++i) state[i] = i+1;
    cout << "state->ket_tensor: " << endl;
    tensor_lattice->VectorizeTensorLattice(false, state);
    tensor_lattice->PrintTensorLattice();
    delete[] state;
   
    cout << endl;
    cout << "5. LeftCanonicalTensorLattice RightCanonicalTensorLattice" << endl;
    cout << "#################################################################" << endl;
    cout << "                    LeftCanonicalTensorLattice                   " << endl;
    cout << "#################################################################" << endl;
    cout << "--------------------------------------------" << endl;
    cout << "Before Left Canonical:" << endl;
    cout << "--------------------------------------------" << endl;
    tensor_lattice->PrintTensorLattice();
    tensor_lattice->LeftCanonicalTensorLattice(-1, 0.0);
    cout << "--------------------------------------------" << endl;
    cout << "After Left Canonical: Ket Tensor" << endl;
    cout << "--------------------------------------------" << endl;
    tensor_lattice->PrintTensorLattice();
    cout << "--------------------------------------------" << endl;
    cout << "After Left Canonical: Canonical Tensor" << endl;
    cout << "--------------------------------------------" << endl;
    tensor_lattice->get_canonical_tensor()->PrintMatrixBlock();

    cout << "#################################################################" << endl;
    cout << "                    RightCanonicalTensorLattice                  " << endl;
    cout << "#################################################################" << endl;
    cout << "--------------------------------------------" << endl;
    cout << "Before Right Canonical:" << endl;
    cout << "--------------------------------------------" << endl;
    tensor_lattice->PrintTensorLattice();
    tensor_lattice->RightCanonicalTensorLattice(5, 0.0);
    cout << "--------------------------------------------" << endl;
    cout << "After Right Canonical: Ket Tensor" << endl;
    cout << "--------------------------------------------" << endl;
    tensor_lattice->PrintTensorLattice();
    cout << "--------------------------------------------" << endl;
    cout << "After Right Canonical: Canonical Tensor" << endl;
    cout << "--------------------------------------------" << endl;
    tensor_lattice->get_canonical_tensor()->PrintMatrixBlock();
    return 0;
}   
