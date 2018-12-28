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
        left_block_[i] = tmp_tensor_lattice->left_block_[i];
        left_dim_[i] = tmp_tensor_lattice->left_dim_[i];
    }

    num_right_block_ = tmp_tensor_lattice->num_right_block_;
    right_block_ = new int[num_right_block_];
    right_dim_ = new int[num_right_block_];
    for(int i=0;i<num_right_block_;++i)
    {
        right_block_[i] = tmp_tensor_lattice->right_block_[i];
        right_dim_[i] = tmp_tensor_lattice->right_dim_[i];
    }

    physics_index_ = new int* [num_left_block_];
    for(int i=0;i<num_left_block_;++i)
    {
        physics_index_[i] = new int[num_right_block_];
        for(int j=0;j<num_right_block_;++j)
        {
            physics_index_[i][j] = tmp_tensor_lattice->physics_index_[i][j];
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
    return num_right_block_;
}

int* RealTensorLattice::get_right_block()
{
    return right_block_;
}

int* RealTensorLattice::get_right_dim()
{
    return right_dim_;
}

int** RealTensorLattice::get_physics_index()
{
    return physics_index_;
}

RealMatrixBlock* RealTensorLattice::get_ket_tensor()
{
    return ket_tensor_;
}

RealMatrixBlock* RealTensorLattice::get_canonical_tensor()
{
    return canonical_tensor_;
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

int RealTensorLattice::ComputePartKetTensorDim(int position)
{
    return ket_tensor_->ComputePartMatrixBlockDim(position);
}

void RealTensorLattice::PrintTensorLattice()
{
    cout << "******************************" << endl;
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
    cout << "ket tensor :" << endl;
    if(ket_tensor_ == nullptr)
        cout << "ket tensor is nullptr!" << endl;
    else
        ket_tensor_->PrintMatrixBlock();
}

void RealTensorLattice::WriteTensorLattice(const char* tensor_lattice_name)
{
    ofstream tensor_lattice_file;

    tensor_lattice_file.open(tensor_lattice_name, ios::binary|ios::out);

    WriteTensorLattice(tensor_lattice_file);

    tensor_lattice_file.close();
}

void RealTensorLattice::WriteTensorLattice(ofstream &tensor_lattice_file)
{
    tensor_lattice_file.write((char*) &physics_dim_, sizeof(int));

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

void RealTensorLattice::ReadTensorLattice(const char* tensor_lattice_name)
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

    tensor_lattice_file.read((char*) &physics_dim_, sizeof(int));

    tensor_lattice_file.read((char*) &num_left_block_, sizeof(int));
    left_block_ = new int[num_left_block_];
    left_dim_ = new int[num_left_block_];
    for(int i=0;i<num_left_block_;++i)
    {
        tensor_lattice_file.read((char*) &left_block_[i], sizeof(int));
        tensor_lattice_file.read((char*) &left_dim_[i], sizeof(int));
    }
    tensor_lattice_file.read((char*) &num_right_block_, sizeof(int));
    right_block_ = new int[num_right_block_];
    right_dim_ = new int[num_right_block_];
    for(int i=0;i<num_right_block_;++i)
    {
        tensor_lattice_file.read((char*) &right_block_[i], sizeof(int));
        tensor_lattice_file.read((char*) &right_dim_[i], sizeof(int));
    }
    physics_index_ = new int* [num_left_block_];
    for(int i=0;i<num_left_block_;++i)
    {
        physics_index_[i] = new int[num_right_block_];
        for(int j=0;j<num_right_block_;j++)
            tensor_lattice_file.read((char*) &physics_index_[i][j], sizeof(int));
    }
    ket_tensor_ = new RealMatrixBlock();
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

void RealTensorLattice::CombineTensorLattice(int leigh, int &num_leigh_block, int* &leigh_block, 
     int* &leigh_dim, RealTensorLattice* expan_tensor_lattice)
{
    int tmp_num_leigh_block[2], *tmp_leigh_block[2], *tmp_leigh_dim[2];
    bool flag;
    if(leigh == 0)
    {
        if(num_right_block_ != expan_tensor_lattice->num_right_block_)
        {
            cout << "num_right_dim_ is not equal to expan_tensor_lattice in CombineTensorLattice!" << endl;
            exit(-1);
        }
        for(int r=0;r<num_right_block_;++r)
        {
            if(right_block_[r] != expan_tensor_lattice->right_block_[r])
            {
                cout << "right_block_ is not equal to expan_tensor_lattice in CimbineTensorLattice!" << endl;
                exit(-1);
            }
        }
    }
    else if(leigh == 1)
    {
        if(num_left_block_ != expan_tensor_lattice->num_left_block_)
        {
            cout << "num_left_dim_ is not equal to expan_tensor_lattice in CombineTensorLattice!" << endl;
            exit(-1);
        }
        for(int l=0;l<num_left_block_;++l)
        {
            if(left_block_[l] != expan_tensor_lattice->left_block_[l])
            {
                cout << "left_block_ is not equal to expan_tensor_lattice in CimbineTensorLattice!" << endl;
                exit(-1);
            }
        }
    }

    if(leigh == 0)
    {
        tmp_num_leigh_block[0] = num_left_block_;
        tmp_leigh_block[0] = left_block_;
        tmp_leigh_dim[0] = left_dim_;
        
        tmp_num_leigh_block[1] = expan_tensor_lattice->num_left_block_;
        tmp_leigh_block[1] = expan_tensor_lattice->left_block_;
        tmp_leigh_dim[1] = expan_tensor_lattice->left_dim_;
    }
    else if(leigh == 1)
    {
        tmp_num_leigh_block[0] = num_right_block_;
        tmp_leigh_block[0] = right_block_;
        tmp_leigh_dim[0] = right_dim_;
        
        tmp_num_leigh_block[1] = expan_tensor_lattice->num_right_block_;
        tmp_leigh_block[1] = expan_tensor_lattice->right_block_;
        tmp_leigh_dim[1] = expan_tensor_lattice->right_dim_;
    }

    num_leigh_block = tmp_num_leigh_block[0];

    // compute new num_leigh_block
    for(int i=0;i<tmp_num_leigh_block[0];++i)
    {
        flag = false;
        for(int j=0;j<tmp_num_leigh_block[1];++j)
        {
            if(tmp_leigh_block[0][i] == tmp_leigh_block[1][j]) 
            {
                flag = true;
                break;
            }
        }
        if(flag == false) num_leigh_block++;
    }
    
    leigh_block = new int[num_leigh_block];
    leigh_dim = new int[num_leigh_block];
    for(int i=0;i<tmp_num_leigh_block[0];++i) leigh_dim[i] = tmp_leigh_dim[0][i];
    // compute new leigh_block and leigh_dim
    int k = tmp_num_leigh_block[0];
    for(int j=0;j<tmp_num_leigh_block[1];++j)
    {
        flag = false;
        for(int i=0;i<tmp_num_leigh_block[0];++i)
        {
            if(tmp_leigh_block[0][i] == tmp_leigh_block[1][j])
            {
                flag = true;
                leigh_dim[i] += tmp_leigh_dim[1][j];
                break;
            }
        }
        if(flag = false)
        {
            leigh_block[k] = tmp_leigh_block[1][j];
            leigh_dim[k] = tmp_leigh_dim[1][j];
            k++;
        }
    }
    // reorder leigh_block ascending depend on leigh_block
    QuickSort(leigh_block, leigh_dim, 0, num_leigh_block-1, 1);
}

void RealTensorLattice::ResetTensorLattice()
{
    ket_tensor_->ResetMatrixBlock();
}

void RealTensorLattice::NormalizeTensorLattice()
{
    ket_tensor_->NormalizeMatrixBlock();
}

void RealTensorLattice::VectorizeTensorLattice(bool direction, double* &state)
{
    ket_tensor_->VectorizeMatrixBlock(direction, state);
}

void RealTensorLattice::ComputeTruncateDim(int max_dim, double canonical_precision, 
     int num_singular_block, int* singular_dim, double** singular_value, int* truncate_dim)
{
    int total_dim = 0;
    double* tmp_singular_value;
    int* tmp_singular_block;
    if(canonical_precision > 1) canonical_precision = 0.0;

    // compute truncate dimension for every block if singular values are 
    // bigger then canonical precisiom
    for(int i=0;i<num_singular_block;++i)
    {
        truncate_dim[i] = 0;
        for(int j=0;j<singular_dim[i];++j)
        {
            if(singular_value[i][j] > canonical_precision) truncate_dim[i]++;
        }
        total_dim += truncate_dim[i];
    }

    // reorder singular values and throw away the smallest one until
    // total dimension is equal to max dimension
    int shift = 0;
    if(total_dim>0 && max_dim>0 && max_dim<total_dim)
    {
        tmp_singular_value = new double[total_dim];
        tmp_singular_block = new int[total_dim];
        for(int i=0;i<num_singular_block;++i)
        {
            for(int j=0;j<truncate_dim[i];++j)
            {
                tmp_singular_value[shift] = singular_value[i][j];
                tmp_singular_block[shift] = i;
                shift++;
            }
        }
        for(int i=0;i<num_singular_block;++i) truncate_dim[i] = 0;
        QuickSort<double>(tmp_singular_value, tmp_singular_block, 0, total_dim-1, 0);
        for(int i=0;i<max_dim;++i) truncate_dim[tmp_singular_block[i]]++;
        delete[] tmp_singular_value;
        delete[] tmp_singular_block;
    }

    // if truncate dimension of one block is equal to zero, change singular
    // value to 0.0 and truncation dimension to 1.
    for(int i=0;i<num_singular_block;++i)
    {
        if(truncate_dim[i] == 0)
        {
            truncate_dim[i] = 1;
            singular_value[i][0] = 0.0;
        }
    }
    
}


void RealTensorLattice::LeftCanonicalTensorLattice(int max_dim, double canonical_precision)
{
    RealMatrix **left_singular_tensor, **right_singular_tensor, *block_tensor,
               *tmp_tensor;
    double **singular_value;
    int *singular_dim, *truncate_dim, block_dim[2], tmp_dim[2], num_same_index, 
        *position_same_index, add_up_left_dim;

    left_singular_tensor = new RealMatrix* [num_right_block_];
    right_singular_tensor = new RealMatrix* [num_right_block_];
    singular_value = new double* [num_right_block_];
    singular_dim = new int[num_right_block_];
    truncate_dim = new int[num_right_block_];
    // do SVD decomposition for every right block number
    //
    for(int r=0;r<num_right_block_;++r)
    {
        ket_tensor_->FindMatrixBlock(num_same_index, position_same_index, 1, r);
        block_dim[0] = 0;
        //block_dim[1] = ket_tensor_->get_matrix_block(position_same_index[0])->get_column();
        block_dim[1] = right_dim_[r];
        for(int p=0;p<num_same_index;++p)
        {
            tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
            block_dim[0] += tmp_tensor->get_row();
        }

        if(block_dim[0] == 0)
        {
            left_singular_tensor[r] = nullptr;
            right_singular_tensor[r] = new RealMatrix(block_dim[1], block_dim[1]);
            singular_dim[r] = 1;
            // new double[0]?
            singular_value[r] = new double[block_dim[1]];
            singular_value[r][0] = 0.0;
        }
        else
        {
            block_tensor = new RealMatrix(block_dim[0], block_dim[1]);
            // put matrices with the same right block number into block_tenosr.
            add_up_left_dim = 0;
            for(int p=0;p<num_same_index;++p)
            {
                tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
                tmp_dim[0] = tmp_tensor->get_row();
                tmp_dim[1] = tmp_tensor->get_column();
                for(int i=0;i<tmp_dim[0];++i) for(int j=0;j<tmp_dim[1];++j)
                {
                    block_tensor->set_matrix_element(add_up_left_dim+i, j, 
                                  tmp_tensor->get_matrix_element(i, j));
                }
                add_up_left_dim += tmp_dim[0];
            }
            // svd
            block_tensor->SVDMatrix(left_singular_tensor[r], right_singular_tensor[r], 
                                    singular_value[r], singular_dim[r]);
            delete block_tensor;
        }

        delete[] position_same_index;
    }
    ComputeTruncateDim(max_dim, canonical_precision, num_right_block_, singular_dim, 
                       singular_value, truncate_dim);

    // left_canonical_tensor => ket_tensor_
    // singular_value*right_canonical_tensor_ => canonical_tensor_
    // match_dim_ : right_dim_ for canonical_tensor_
    delete[] match_dim_;
    delete canonical_tensor_;
    
    match_dim_ = new int[num_right_block_];
    canonical_tensor_ = new RealMatrixBlock(num_right_block_);
    for(int r=0;r<num_right_block_;++r)
    {
        if(left_singular_tensor[r] != nullptr)
        {
            ket_tensor_->FindMatrixBlock(num_same_index, position_same_index, 1, r);
            add_up_left_dim = 0;

            for(int p=0;p<num_same_index;++p)
            {
                // tmp_tensor is a pointer which points to a matrix from ket_tensor
                tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
                tmp_tensor->ChangeMatrix(1, truncate_dim[r]);

                tmp_dim[0] = tmp_tensor->get_row();
                tmp_dim[1] = tmp_tensor->get_column();
                for(int i=0;i<tmp_dim[0];++i) for(int j=0;j<tmp_dim[1];++j)
                {
                    tmp_tensor->set_matrix_element(i, j, 
                    left_singular_tensor[r]->get_matrix_element(add_up_left_dim+i, j));
                }
                add_up_left_dim += tmp_dim[0];
            }
            delete[] position_same_index;
        }

        right_dim_[r] = truncate_dim[r];
        block_dim[1] = right_singular_tensor[r]->get_column();
        match_dim_[r] = block_dim[1];
        block_tensor = new RealMatrix(truncate_dim[r], block_dim[1]);
        for(int i=0;i<truncate_dim[r];++i) for(int j=0;j<block_dim[1];++j)
            block_tensor->set_matrix_element(i, j, singular_value[r][i]*
                    right_singular_tensor[r]->get_matrix_element(i, j));
        
        canonical_tensor_->set_matrix_block(r, block_tensor);
        
        delete left_singular_tensor[r];
        delete right_singular_tensor[r];
        delete[] singular_value[r];
        delete block_tensor;
    }
    delete[] left_singular_tensor;
    delete[] right_singular_tensor;
    delete[] singular_value;
    delete[] singular_dim;
    delete[] truncate_dim;
}

void RealTensorLattice::LeftCanonicalTensorLattice(int* &singular_dim, double** &singular_value)
{
    RealMatrix **left_singular_tensor, **right_singular_tensor, *block_tensor,
               *tmp_tensor;
    int *truncate_dim, block_dim[2], tmp_dim[2], num_same_index, 
        *position_same_index, add_up_left_dim;

    left_singular_tensor = new RealMatrix* [num_right_block_];
    right_singular_tensor = new RealMatrix* [num_right_block_];
    singular_value = new double* [num_right_block_];
    singular_dim = new int[num_right_block_];
    truncate_dim = new int[num_right_block_];
    // do SVD decomposition for every right block number
    //
    for(int r=0;r<num_right_block_;++r)
    {
        ket_tensor_->FindMatrixBlock(num_same_index, position_same_index, 1, r);
        block_dim[0] = 0;
        //block_dim[1] = ket_tensor_->get_matrix_block(position_same_index[0])->get_column();
        block_dim[1] = right_dim_[r];
        for(int p=0;p<num_same_index;++p)
        {
            tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
            block_dim[0] += tmp_tensor->get_row();
        }

        if(block_dim[0] == 0)
        {
            left_singular_tensor[r] = nullptr;
            right_singular_tensor[r] = new RealMatrix(block_dim[1], block_dim[1]);
            singular_dim[r] = 1;
            // new double[0]?
            singular_value[r] = new double[block_dim[1]];
            singular_value[r][0] = 0.0;
        }
        else
        {
            block_tensor = new RealMatrix(block_dim[0], block_dim[1]);
            // put matrices with the same right block number into block_tenosr.
            add_up_left_dim = 0;
            for(int p=0;p<num_same_index;++p)
            {
                tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
                tmp_dim[0] = tmp_tensor->get_row();
                tmp_dim[1] = tmp_tensor->get_column();
                for(int i=0;i<tmp_dim[0];++i) for(int j=0;j<tmp_dim[1];++j)
                {
                    block_tensor->set_matrix_element(add_up_left_dim+i, j, 
                                  tmp_tensor->get_matrix_element(i, j));
                }
                add_up_left_dim += tmp_dim[0];
            }
            // svd
            block_tensor->SVDMatrix(left_singular_tensor[r], right_singular_tensor[r], 
                                    singular_value[r], singular_dim[r]);
            delete block_tensor;
        }

        delete[] position_same_index;
    }
    ComputeTruncateDim(-1, 0.0, num_right_block_, singular_dim, 
                       singular_value, truncate_dim);

    // left_canonical_tensor => ket_tensor_
    // singular_value*right_canonical_tensor_ => canonical_tensor_
    // match_dim_ : right_dim_ for canonical_tensor_
    delete[] match_dim_;
    delete canonical_tensor_;
    
    match_dim_ = new int[num_right_block_];
    canonical_tensor_ = new RealMatrixBlock(num_right_block_);
    for(int r=0;r<num_right_block_;++r)
    {
        if(left_singular_tensor[r] != nullptr)
        {
            ket_tensor_->FindMatrixBlock(num_same_index, position_same_index, 1, r);
            add_up_left_dim = 0;

            for(int p=0;p<num_same_index;++p)
            {
                // tmp_tensor is a pointer which points to a matrix from ket_tensor
                tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
                tmp_tensor->ChangeMatrix(1, truncate_dim[r]);

                tmp_dim[0] = tmp_tensor->get_row();
                tmp_dim[1] = tmp_tensor->get_column();
                for(int i=0;i<tmp_dim[0];++i) for(int j=0;j<tmp_dim[1];++j)
                {
                    tmp_tensor->set_matrix_element(i, j, 
                    left_singular_tensor[r]->get_matrix_element(add_up_left_dim+i, j));
                }
                add_up_left_dim += tmp_dim[0];
            }
            delete[] position_same_index;
        }

        right_dim_[r] = truncate_dim[r];
        block_dim[1] = right_singular_tensor[r]->get_column();
        match_dim_[r] = block_dim[1];
        block_tensor = new RealMatrix(truncate_dim[r], block_dim[1]);
        for(int i=0;i<truncate_dim[r];++i) for(int j=0;j<block_dim[1];++j)
            block_tensor->set_matrix_element(i, j, singular_value[r][i]*
                    right_singular_tensor[r]->get_matrix_element(i, j));
        
        canonical_tensor_->set_matrix_block(r, block_tensor);
        
        delete left_singular_tensor[r];
        delete right_singular_tensor[r];
        delete block_tensor;
    }
    delete[] left_singular_tensor;
    delete[] right_singular_tensor;
    delete[] truncate_dim;
}


void RealTensorLattice::RightCanonicalTensorLattice(int max_dim, double canonical_precision)
{
    RealMatrix **left_singular_tensor, **right_singular_tensor, *block_tensor,
               *tmp_tensor;
    double **singular_value;
    int *singular_dim, *truncate_dim, block_dim[2], tmp_dim[2], num_same_index, 
        *position_same_index, add_up_right_dim;

    left_singular_tensor = new RealMatrix* [num_left_block_];
    right_singular_tensor = new RealMatrix* [num_left_block_];
    singular_value = new double* [num_left_block_];
    singular_dim = new int[num_left_block_];
    truncate_dim = new int[num_left_block_];
    // do SVD decomposition for every left block number
    //
    for(int l=0;l<num_left_block_;++l)
    {
        ket_tensor_->FindMatrixBlock(num_same_index, position_same_index, 0, l);
        block_dim[0] = left_dim_[l];
        block_dim[1] = 0;
        for(int p=0;p<num_same_index;++p)
        {
            tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
            block_dim[1] += tmp_tensor->get_column();
        }

        if(block_dim[1] == 0)
        {
            left_singular_tensor[l] = new RealMatrix(block_dim[0], block_dim[0]);
            right_singular_tensor[l] = nullptr;
            singular_dim[l] = 1;
            // new double[0]?
            singular_value[l] = new double[block_dim[0]];
            singular_value[l][0] = 0.0;
        }
        else
        {
            block_tensor = new RealMatrix(block_dim[0], block_dim[1]);
            // put matrices with the same left block number into block_tenosr.
            add_up_right_dim = 0;
            for(int p=0;p<num_same_index;++p)
            {
                tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
                tmp_dim[0] = tmp_tensor->get_row();
                tmp_dim[1] = tmp_tensor->get_column();
                for(int i=0;i<tmp_dim[0];++i) for(int j=0;j<tmp_dim[1];++j)
                {
                    block_tensor->set_matrix_element(i, add_up_right_dim+j, 
                                  tmp_tensor->get_matrix_element(i, j));
                }
                add_up_right_dim += tmp_dim[1];
            }
            // svd
            block_tensor->SVDMatrix(left_singular_tensor[l], right_singular_tensor[l], 
                                    singular_value[l], singular_dim[l]);
            delete block_tensor;
        }

        delete[] position_same_index;
    }
    ComputeTruncateDim(max_dim, canonical_precision, num_right_block_, singular_dim, 
                       singular_value, truncate_dim);
    // right_canonical_tensor => ket_tensor_
    // left_canonical_tensor_*singular_value => canonical_tensor_
    // match_dim_ : left_dim_ for canonical_tensor_
    delete[] match_dim_;
    delete canonical_tensor_;
    
    match_dim_ = new int[num_left_block_];
    canonical_tensor_ = new RealMatrixBlock(num_left_block_);
    for(int l=0;l<num_left_block_;++l)
    {
        if(right_singular_tensor[l] != nullptr)
        {
            ket_tensor_->FindMatrixBlock(num_same_index, position_same_index, 0, l);
            add_up_right_dim = 0;

            for(int p=0;p<num_same_index;++p)
            {
                // tmp_tensor is a pointer which points to a matrix from ket_tensor
                tmp_tensor = ket_tensor_->get_matrix_block(position_same_index[p]);
                tmp_tensor->ChangeMatrix(0, truncate_dim[l]);

                tmp_dim[0] = tmp_tensor->get_row();
                tmp_dim[1] = tmp_tensor->get_column();
                for(int i=0;i<tmp_dim[0];++i) for(int j=0;j<tmp_dim[1];++j)
                {
                    tmp_tensor->set_matrix_element(i, j, 
                    right_singular_tensor[l]->get_matrix_element(i, add_up_right_dim+j));
                }
                add_up_right_dim += tmp_dim[1];
            }
            delete[] position_same_index;
        }

        left_dim_[l] = truncate_dim[l];
        block_dim[0] = left_singular_tensor[l]->get_row();
        match_dim_[l] = block_dim[0];
        block_tensor = new RealMatrix(block_dim[0], truncate_dim[l]);
        for(int i=0;i<block_dim[0];++i) for(int j=0;j<truncate_dim[l];++j)
            block_tensor->set_matrix_element(i, j, left_singular_tensor[l]->get_matrix_element(i, j)
                    *singular_value[l][j]);
        
        canonical_tensor_->set_matrix_block(l, block_tensor);
        
        delete left_singular_tensor[l];
        delete right_singular_tensor[l];
        delete[] singular_value[l];
        delete block_tensor;
    }
    delete[] left_singular_tensor;
    delete[] right_singular_tensor;
    delete[] singular_value;
    delete[] singular_dim;
    delete[] truncate_dim;
}


void RealTensorLattice::LeftMergeTensorLattice(RealTensorLattice* tmp_tensor_lattice)
{
    RealMatrix *merge_tensor, *tmp_tensor[2];
    int num_same_index, *position_same_index;

    if(tmp_tensor_lattice->get_canonical_tensor() == nullptr)
    {
        cout << "Canonical tensor is nullptr in LeftMergeTensorLattice!";
        exit(-1);
    }
    for(int l=0;l<num_left_block_;++l)
    {
        ket_tensor_->FindMatrixBlock(num_same_index, position_same_index, 0, l);
        // direct get protected members(why can't)
        tmp_tensor[0] = tmp_tensor_lattice->canonical_tensor_->get_matrix_block(l);
        for(int p=0;p<num_same_index;++p)
        {
            tmp_tensor[1] = ket_tensor_->get_matrix_block(position_same_index[p]);
            merge_tensor = tmp_tensor[0]->MultiplyToMatrix(tmp_tensor[1]);
            tmp_tensor[1]->ReplaceMatrix(merge_tensor);
            delete merge_tensor;
        }
        left_dim_[l] =tmp_tensor[0]->get_row();
        delete[] position_same_index;
    }
    NormalizeTensorLattice();
    delete[] tmp_tensor_lattice->match_dim_;
    delete tmp_tensor_lattice->canonical_tensor_;
    tmp_tensor_lattice->match_dim_ = nullptr;
    tmp_tensor_lattice->canonical_tensor_ = nullptr;
}


void RealTensorLattice::RightMergeTensorLattice(RealTensorLattice* tmp_tensor_lattice)
{
    RealMatrix *merge_tensor, *tmp_tensor[2];
    int num_same_index, *position_same_index;

    if(tmp_tensor_lattice->get_canonical_tensor() == nullptr)
    {
        cout << "Canonical tensor is nullptr in RightMergeTensorLattice!";
        exit(-1);
    }
    for(int r=0;r<num_right_block_;++r)
    {
        ket_tensor_->FindMatrixBlock(num_same_index, position_same_index, 1, r);
        // direct get protected members
        tmp_tensor[1] = tmp_tensor_lattice->canonical_tensor_->get_matrix_block(r);
        for(int p=0;p<num_same_index;++p)
        {
            tmp_tensor[0] = ket_tensor_->get_matrix_block(position_same_index[p]);
            merge_tensor = tmp_tensor[0]->MultiplyToMatrix(tmp_tensor[1]);
            tmp_tensor[0]->ReplaceMatrix(merge_tensor);
            delete merge_tensor;
        }
        right_dim_[r] =tmp_tensor[1]->get_column();
        delete[] position_same_index;
    }
    NormalizeTensorLattice();
    delete[] tmp_tensor_lattice->match_dim_;
    delete tmp_tensor_lattice->canonical_tensor_;
    tmp_tensor_lattice->match_dim_ = nullptr;
    tmp_tensor_lattice->canonical_tensor_ = nullptr;
}

void RealTensorLattice::LeftExpanTensorLattice(RealTensorLattice* tmp_tensor_lattice, 
        RealTensorLattice* expan_tensor_lattice)
{
    RealMatrix *origin_tensor, *tmp_tensor, *expan_tensor;
    int position, enable_expan[2], position_expan[2];

    for(int r=0;r<num_right_block_;++r)
    {
        enable_expan[0] = -1;
        for(int i=0;i<tmp_tensor_lattice->num_right_block_;++i)
        {
            if(right_block_[i] == tmp_tensor_lattice->right_block_[i])
            {
                enable_expan[0] = i;
                break;
            }
        }

        enable_expan[1] = -1;
        if(expan_tensor_lattice != nullptr)
        {
            for(int i=0;i<expan_tensor_lattice->num_right_block_;++i)
            {
                if(right_block_[i] == expan_tensor_lattice->right_block_[i])
                {
                    enable_expan[1] = i;
                    break;
                }
            }
        }

        for(int l=0;l<num_left_block_;++l)
        {
            position = ket_tensor_->FindMatrixBlock(l, r);
            if(position != -1)
            {
                origin_tensor = ket_tensor_->get_matrix_block(position);

                if(enable_expan[0] != -1)
                {
                    position_expan[0] = tmp_tensor_lattice->ket_tensor_->FindMatrixBlock(l, enable_expan[0]);
                    if(position_expan[0] != -1)
                    {
                        tmp_tensor = tmp_tensor_lattice->ket_tensor_->get_matrix_block(position_expan[0]);
                        origin_tensor->AddToMatrix(tmp_tensor);
                    }
                }

                if(enable_expan[1] != -1)
                {
                    position_expan[1] = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(l, enable_expan[1]);
                    if(position_expan[0] != -1)
                    {
                        expan_tensor = expan_tensor_lattice->ket_tensor_->get_matrix_block(position_expan[1]);
                        origin_tensor->ExpanMatrix(1, expan_tensor);
                    }

                }

                if(enable_expan[0]!=-1 && enable_expan[1]==-1)
                {
                    cout << "enable_expan[0]!=-1 && enable_expan[1]==-1 in LeftExpanTensorLattice!" << endl;
                }
                if(enable_expan[0]==-1 && enable_expan[1]==-1)
                {
                    cout << "enable_expan[0]==-1 && enable_expan[1]==-1 in LeftExpanTensorLattice!" << endl;
                }
            }
        }
    }
}


void RealTensorLattice::RightExpanTensorLattice(RealTensorLattice* tmp_tensor_lattice, 
        RealTensorLattice* expan_tensor_lattice)
{
    RealMatrix *origin_tensor, *tmp_tensor, *expan_tensor;
    int position, enable_expan[2], position_expan[2];

    for(int l=0;l<num_left_block_;++l)
    {
        enable_expan[0] = -1;
        for(int i=0;i<tmp_tensor_lattice->num_left_block_;++i)
        {
            if(left_block_[i] == tmp_tensor_lattice->left_block_[i])
            {
                enable_expan[0] = i;
                break;
            }
        }

        enable_expan[1] = -1;
        if(expan_tensor_lattice != nullptr)
        {
            for(int i=0;i<expan_tensor_lattice->num_left_block_;++i)
            {
                if(left_block_[i] == expan_tensor_lattice->left_block_[i])
                {
                    enable_expan[1] = i;
                    break;
                }
            }
        }

        for(int r=0;r<num_right_block_;++r)
        {
            position = ket_tensor_->FindMatrixBlock(l, r);
            if(position != -1)
            {
                origin_tensor = ket_tensor_->get_matrix_block(position);

                if(enable_expan[0] != -1)
                {
                    position_expan[0] = tmp_tensor_lattice->ket_tensor_->FindMatrixBlock(enable_expan[0], r);
                    if(position_expan[0] != -1)
                    {
                        tmp_tensor = tmp_tensor_lattice->ket_tensor_->get_matrix_block(position_expan[0]);
                        origin_tensor->AddToMatrix(tmp_tensor);
                    }
                }

                if(enable_expan[1] != -1)
                {
                    position_expan[1] = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(enable_expan[1], r);
                    if(position_expan[0] != -1)
                    {
                        expan_tensor = expan_tensor_lattice->ket_tensor_->get_matrix_block(position_expan[1]);
                        origin_tensor->ExpanMatrix(0, expan_tensor);
                    }

                }

                if(enable_expan[0]!=-1 && enable_expan[1]==-1)
                {
                    cout << "enable_expan[0]!=-1 && enable_expan[1]==-1 in RightExpanTensorLattice!" << endl;
                }
                if(enable_expan[0]==-1 && enable_expan[1]==-1)
                {
                    cout << "enable_expan[0]==-1 && enable_expan[1]==-1 in RightExpanTensorLattice!" << endl;
                }
            }
        }
    }
}



