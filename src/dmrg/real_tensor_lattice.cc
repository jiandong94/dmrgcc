#include "dmrg/real_tensor_lattice.h"


RealTensorLattice::RealTensorLattice()
{
    physics_dim_ = 0;

    num_left_block_ = 0;
    left_block_ = nullptr;
    left_dim_ = nullptr;

    num_right_block_ = 0;
    right_block_ = nullptr;
    right_dim_ = nullptr;

    physics_index_ = nullptr;

    ket_tensor_ = nullptr;

    match_dim_ =  nullptr;
    canonical_tensor_ = nullptr;
}


RealTensorLattice::RealTensorLattice(int physics_dim)
{
    physics_dim_ = physics_dim;

    num_left_block_ = 0;
    left_block_ = nullptr;
    left_dim_ = nullptr;

    num_right_block_ = 0;
    right_block_ = nullptr;
    right_dim_ = nullptr;

    physics_index_ = nullptr;

    ket_tensor_ = nullptr;

    match_dim_ =  nullptr;
    canonical_tensor_ = nullptr;
}


RealTensorLattice::RealTensorLattice(RealTensorLattice* tmp_tensor_lattice)
{
    physics_dim_ = tmp_tensor_lattice->physics_dim_;

    num_left_block_ = tmp_tensor_lattice->num_left_block_;
    left_block_ = new int[num_left_block_];
    left_dim_ = new int[num_left_block_];
    for(int i=0;i<num_left_block_;++i)
    {
        left_block_ = tmp_tensor_lattice->left_block_[i];
        left_dim_ = tmp_tensor_lattice->left_dim_[i];
    }

    num_right_block_ = tmp_tensor_lattice->num_right_block_;
    right_block_ = new int[num_right_block_];
    right_dim_ = new int[num_right_block_];
    for(int i=0;i<num_right_block_;++i)
    {
        right_block_ = tmp_tensor_lattice->right_block_[i];
        right_dim_ = tmp_tensor_lattice->right_dim_[i];
    }

    physics_index_ = new int* [num_left_block_];
    for(int i=0;i<num_left_block_;++i)
    {
        physics_index_[i] = new int[num_right_block_];
        for(int j=0;j<num_right_block_;++j)
        {
            phycics_index_[i][j] = tmp_tensor_lattice->physics_index_[i][j];
        }
    }
    

    ket_tensor_ = new RealMatrixBlock(tmp_tensor_lattice->ket_tensor_);

    match_dim_ =  nullptr;
    canonical_tensor_ = nullptr;
}

RealTensorLattice::~RealTensorLattice()
{
    delete[] left_block_;
    delete[] left_dim_;

    delete[] right_block_;
    delete[] right_dim_;

    for(int i=0;i<num_left_block_;++i) delete[] physics_index_[i];
    delete[] physics_index_;

    delete ket_tensor_;

    delete[] match_dim_;
    delete canonical_tensor_;
}

int RealTensorLattice::get_physics_dim()
{
    return physics_dim_;
}

int RealTensorLattice::get_num_left_block()
{
    return num_left_block_;
}

int* RealTensorLattice::get_left_block()
{
    return left_block_;
}

int* RealTensorLattice::get_left_dim()
{
    return left_dim_;
}

int RealTensorLattice::get_num_right_block()
{
    return num_left_block_;
}

int* RealTensorLattice::get_right_block()
{
    return right_block_;
}

int* RealTensorLattice::get_right_dim_()
{
    return right_dim_;
}

int** RealTensorLattice::get_physics_index()
{
    return physics_index_;
}

int RealTensorLattice::ComputeLatticeDim(int leigh)
{
    int bond_dim = 0;
    if(leigh == 0)
    {
        for(int i=0;i<num_left_block_;++i) bond_dim+=left_dim_[i];
    }
    else if(leigh == 1)
    {
        for(int i=0;i<num_right_block_;++i) bond_dim+=right_dim_[i];
    }
    return bond_dim;
}

int RealTensorLattice::ComputeKetTensorDim()
{
    return ket_tensor_->ComputeMatrixBlockDim();
}

int RealTensorLattice::ComputePartKetTensorDim()
{
    return ket_tensor_->ComputePartMatrixBlockDim();
}

void RealTensorLattice::PrintTensorLattice()
{
    cout << "******************************" << endl;
    cout << "Print Tensor Lattice" << endl;
    cout << "number of left  block = " << num_left_block_ << ",";
    cout << "number of right block = " << num_right_block_ << endl;

    cout << "left  block: ";
    for(int i=0;i<num_left_block_;++i) cout << left_block_[i] << ", ";
    cout << endl;
    cout << "left  dim  : ";
    for(int i=0;i<num_left_block_;++i) cout << left_dim_[i] << ", ";
    cout << endl;

    cout << "right block: ";
    for(int i=0;i<num_right_block_;++i) cout << right_block_[i] << ", ";
    cout << endl;
    cout << "right dim  : ";
    for(int i=0;i<num_right_block_;++i) cout << right_dim_[i] << ", ";
    cout << endl;

    ket_tensor_->PrintMatrixBlock();
}

void RealTensorLattice::WriteTensorLattice(char* tensor_lattice_name)
{
    ofstream tensor_lattice_file;

    tensor_lattice_file.open(tensor_lattice_name, ios::binary|ios::out);

    WriteTensorLattice(tensor_lattice_file);

    tensor_lattice_file.close();
}

void RealTensorLattice::WriteTensorLattice(ofstream &tensor_lattice_file)
{
    tensor_lattice_file.write((char*) &physice_dim_, sizeof(int));

    tensor_lattice_file.write((char*) &num_left_block_, sizeof(int));
    for(int i=0;i<num_left_block_;++i)
    {
        tensor_lattice_file.write((char*) &left_block_[i], sizeof(int));
        tensor_lattice_file.write((char*) &left_dim_[i], sizeof(int));
    }
    tensor_lattice_file.write((char*) &num_right_block_, sizeof(int));
    for(int i=0;i<num_right_block_;++i)
    {
        tensor_lattice_file.write((char*) &right_block_[i], sizeof(int));
        tensor_lattice_file.write((char*) &right_dim_[i], sizeof(int));
    }

    for(int i=0;i<num_left_block_;++i) for(int j=0;j<num_right_block_;++j)
        tensor_lattice_file.write((char*) &physics_index_[i][j], sizeof(int));

    ket_tensor_->WriteMatrixBlock(tensor_lattice_file);
}

void RealTensorLattice::ReadTensorLattice(char* tensor_lattice_name)
{
    ifstream tensor_lattice_file;

    tensor_lattice_file.open(tensor_lattice_name, ios::binary|ios::in);

    ReadTensorLattice(tensor_lattice_file);

    tensor_lattice_file.close();
}

void RealTensorLattice::ReadTensorLattice(ifstream &tensor_lattice_file)
{
    delete[] left_block_;
    delete[] left_dim_;

    delete[] right_block_;
    delete[] right_dim_;

    for(int i=0;i<num_left_block_;++i) delete[] physics_index_[i];
    delete[] physics_index_;

    delete ket_tensor_;

    delete[] match_dim_;
    delete canonical_tensor_;

    tensor_lattice_file.read((char*) &physice_dim_, sizeof(int));

    tensor_lattice_file.read((char*) &num_left_block_, sizeof(int));
    left_block_ = new int[num_left_block_];
    left_dim_ = new int[num_left_dim_];
    for(int i=0;i<num_left_block_;++i)
    {
        tensor_lattice_file.read((char*) &left_block_[i], sizeof(int));
        tensor_lattice_file.read((char*) &left_dim_[i], sizeof(int));
    }
    tensor_lattice_file.read((char*) &num_right_block_, sizeof(int));
    right_block_ = new int[num_right_block_];
    right_dim_ = new int[num_right_dim_];
    for(int i=0;i<num_right_block_;++i)
    {
        tensor_lattice_file.read((char*) &right_block_[i], sizeof(int));
        tensor_lattice_file.read((char*) &right_dim_[i], sizeof(int));
    }

    physics_index_ = new int* [num_left_block_];
    for(int i=0;i<num_left_block_;++i)
    {
        physics_index_[i] = new int[num_right_block_];
        tensor_lattice_file.read((char*) &physics_index_[i][j], sizeof(int));
    }
    ket_tensor_->ReadMatrixBlock(tensor_lattice_file);
}

void RealTensorLattice::DefineTensorLattice(int num_left_block, int num_right_block, 
     int* left_block, int* right_block, int* left_dim, int* right_dim)
{
    delete[] left_block_;
    delete[] left_dim_;

    delete[] right_block_;
    delete[] right_dim_;

    for(int i=0;i<num_left_block_;++i) delete[] physics_index_[i];
    delete[] physics_index_;

    num_left_block_ = num_left_block;
    left_block_ = new int[num_left_block_];
    left_dim_ = new int[num_left_block_];
    for(int i=0;i<num_left_block_;++i)
    {
        left_block_[i] = left_block[i];
        left_dim_[i] = left_dim[i];
    }
    
    num_right_block_ = num_right_block;
    right_block_ = new int[num_right_block_];
    right_dim_ = new int[num_right_block_];
    for(int i=0;i<num_right_block_;++i)
    {
        right_block_[i] = right_block[i];
        right_dim_[i] = right_dim[i];
    }

    physics_index_ = new int* [num_left_block_];
    for(int i=0;i<num_left_block_;++i)
    {
        physics_index_[i] = new int[num_right_block_];
        for(int j=0;j<num_right_block_;++j)
        {
            physics_index_[i][j] = -1;
        }
    }
}


void RealTensorLattice::DefineTensorLattice(int num_block, int* left_index, int* right_index, 
     int* physics_index, int row, int column)
{
    RealMatrix* tmp_matrix;
    
    // matrices are deleted?
    delete ket_tensor_;

    ket_tensor_ = new RealMatrixBlock(num_block, left_index, right_index);
    tmp_matrix = new RealMatrix(row, column);
    for(int i=0;i<num_block;++i)
    {
        physics_index_[left_index[i]][right_index[i]] = physics_index[i];
        tmp_matrix->RandomMatrix();
        ket_tensor_->set_matrix_block(i, tmp_matrix);
    }
    delete tmp_matrix;
}

void RealTensorLattice::CombineTensorLattice()


void RealTensorLattice::ResetTensorLattice();
{
    ket_tensor_->ResetMatrixBlock();
}

void RealTensorLattice::NormalizeTensorLattice()
{

}



















