#include "dmrg/real_tensor_contraction.h"

RealTensorContraction::RealTensorContraction(int left_bond, int right_bond)
{
    left_bond_ = left_bond;
    right_bond_ = right_bond;

    left_contraction_tensor_ = RealMatrixBlock* [left_bond_];
    for(int i=0;i<left_bond_;++i)
        left_contraction_tensor_[i] = nullptr;

    right_contraction_tensor_ = RealMatrixBlock* [right_bond_];
    for(int i=0;i<right_bond_;++i)
        right_contraction_tensor_[i] = nullptr;
}

RealTensroContraction::~RealTensorContraction()
{
    for(int i=0;i<left_bond_;++i)
        delete left_contraction_tensor_[i];
    delete[] left_contraction_tensor_;
    
    for(int i=0;i<right_bond_;++i)
        delete right_contraction_tensor_[i];
    delete[] right_contraction_tensor_;
}

int RealTensorContraction::get_left_bond_()
{
    return left_bond_;
}

int RealTensorContraction::get_right_bond_()
{
    return right_bond_;
}

RealMatrixBlock** RealTensorContraction::get_left_contraction_tensor()
{
    return left_contraction_tensor_;
}

RealMatrixBlock** RealTensorContraction::get_right_contraction_tensor()
{
    return right_contraction_tensor_;
}

void RealTensorContraction::PrintTensorContraction()
{
    cout << "left bond = " << left_bond_ << endl;
    cout << "right bond = " << right_bond_ << endl;
    
    cout << "left contraction tensor L : " << endl;
    for(int i=0;i<left_bond_;++i)
    {
        cout << "left bond index = " << i << endl;
        left_contraction_tensor_[i]->PrintMatrixBlock();
    }
    
    cout << "right contraction tensor R : " << endl;
    for(int i=0;i<right_bond_;++i)
    {
        cout << "right bond index = " << i << endl;
        right_contraction_tensor_[i]->PrintMatrixBlock();
    }
}

void RealTensorContraction::WriteTensorContraction(int leigh, char* tensor_contraction_name)
{
    ofstream tensor_contraction_file;
    tensor_contraction_file.open(tensor_space_name, ios::binary|ios::out);

    if(leigh == 0)
    {
        for(int i=0;i<left_bond_;++i)
            left_contraction_tensor_[i]->WriteMatrixBlock(tensor_contraction_file);
    }
    else if(leigh == 1)
    {
        for(int i=0;i<right_bond_;++i)
            right_contraction_tensor_[i]->WriteMatrixBlock(tensor_contraction_file);
    }

    tensor_contraction_file.close();
}

void RealTensorContraction::ReadTensorContraction(int leigh, char* tensor_contraction_name)
{
    ifstream tensor_contraction_file;
    tensor_contraction_file.open(tensor_space_name, ios::binary|ios::in);

    if(leigh == 0)
    {
        for(int i=0;i<left_bond_;++i)
            left_contraction_tensor_[i]->ReadMatrixBlock(tensor_contraction_file);
    }
    else if(leigh == 1)
    {
        for(int i=0;i<right_bond_;++i)
            right_contraction_tensor_[i]->ReadMatrixBlock(tensor_contraction_file);
    }

    tensor_contraction_file.close();
}

void RealTensorContraction::DefineTensorContraction(int leigh)
{
    RealMatrix tmp_tensor;
    
    tmp_tensor = new RealMatrix(1,1);
    tmp_tensor->set_matrix_element(0,0,1);

    if(leigh == 0)
    {
        for(int i=0;i<left_bond_;++i)
        {
            if(left_contraction_tensor_[i] != nullptr) delete left_contraction_tensor_[i];
            left_contraction_tensor_[i] = new RealMatrixBlock(1);
            left_contraction_tensor_[i]->set_matrix_block(0, tmp_tensor);
        }
    }
    else if(leigh == 1)
    {
        for(int i=0;i<right_bond_;++i)
        {
            if(right_contraction_tensor_[i] != nullptr) delete right_contraction_tensor_[i];
            right_contraction_tensor_[i] = new RealMatrixBlock(1);
            right_contraction_tensor_[i]->set_matrix_block(0, tmp_tensor);
        }
    }
    delete tmp_tensor;
}

void RealTensorContraction::ResetTensorContraction(int leigh)
{
    if(leigh == 0)
    {
        for(int i=0;i<left_bond_;++i)
            left_contraction_tensor_[i]->ResetMatrixBlock();
    }
    else if(leigh == 1)
    {
        for(int i=0;i<right_bond_;++i)
            right_contraction_tensor_[i]->ResetMatrixBlock();
    }
}

/*
//L^{o1, k, l}_{b_0, b_1}          L^{o2, i, j}_{a_0, a_1}
//                                          
// ______                           ______ --- ___  
//|      idx:k, dim:b_0            |      - B -   idx:i, dim:a_0
//|                                |       ---  
//|                                |        |p1
//|______                  _____\  |______-----___
//|      idx:o_1                /  |      - W -   idx:o_2
//|                                |      -----            
//|                                |        |p2
//|______                          |______ --- ___
//       idx:l, dim:b_1                   - K -   idx:j, dim:a_1
//                                         ---  
*/                                       
void RealTensorContraction::LeftComputeTensorContraction(RealTensorLattice* tensor_lattice, 
        RealTensorOperator* tensor_operator, RealTensorContraction* tensor_contraction)
{
    RealMatrix *first_tensor, *second_tensor, *tmp_tensor[3], *bra_tensor;
    bool ***zero_flag;
    int *num_same_index, **position_same_index, num_block, *left_index, *right_index;
    int left_bond, right_bond, num_left_block, num_right_block;
    int p1, p2, idx, jdx, kdx, ldx, index;
    double element;

    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();
    
    zero_flag = new int** [right_bond];
    for(int o=0;o<right_bond;++o)
    {
        zero_flag[o] = new int* [num_right_block];
        for(int i=0;i<num_right_block;++i)
        {
            zero_flag[o][i] = new int[num_right_block];
            for(int j=0;j<num_right_block;++j)
                zero_flag[o][i][j] = false;
        }
    }

    // for each left index(block), find all possible right indiced(blocks)
    num_same_index = new int[num_left_block];
    position_same_index = new int* [num_left_block];
    for(int i=0;i<num_left_block;++i)
        tensor_lattice->ket_tensor_->FindMatrixBlock(num_same_index[i], position_same_index[i], 0, i);

    // if zero_flag[o2][i][j] is true, 
    // then left_contraction_tensor[o]->find_matrix_block(i, j) is nonezero
    for(int o1=0;o1<left_bond;++o1)
    for(int i=0;i<left_contraction_tensor_[o1]->get_num_block();++i)
    {
        kdx = left_contracted_tensor_[o1]->left_index_[i];
        ldx = left_contracted_tensor_[o1]->right_index_[i];
        tmp_tensor[0] = left_contracted_tensor_[o1]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int k=0;k<num_same_index[kdx];++k)
        {
            idx = tensor_lattice->ket_tensor_->right_index_[position_same_index[kdx][k]];
            p1 = tensor_lattice->ket_tensor_->physics_index_[kdx][idx];
            if(tensor_operator->LeftCheckZero(o1, p1) == false) continue; //zero
            for(int l=0;l<num_same_index[ldx];++l)
            {
                jdx = tensor_lattice->ket_tensor_->right_index_[position_same_index[ldx][j]];
                p2 = tensor_lattice->ket_tensor_->physics_index_[ldx][jdx];
                for(int o2=0;o2<right_bond;++o2)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_element(p1, p2);
                    if(element != 0.0) zero_flag[o2][idx][jdx] = true;
                }
            }
        }
    }

    // generate left_contraction_tensor accronding to zero_flag,
    // (num_block, left_index, right_index)
    for(int o2=0;o2<right_bond;++o2)
    {
        num_block = 0;
        for(int i=0;i<num_right_block;++i) for(int j=0;j<num_right_block;++j)
        {
            if(zero_flag[o2][i][j] == true)
            {
                num_block++;
            }
        }

        if(num_block == 0)
        {
            left_index = nullptr;
            right_inex = nullptr;
        }
        else
        {
            index = 0;
            left_index = new int[num_block];
            right_index = new int[num_block];
            for(int i=0;i<num_right_block;++i) for(int j=0;j<num_right_block;++j)
            {
                if(zero_flag[o2][i][j] == true) 
                {
                    left_index[index] = i;
                    right_index[index] = j;
                    index++;
                }
            }
        }

        // if right_bond > left_bond?
        delete tensor_contraction->left_contraction_tensor[o2];
        tensor_contraction->left_contraction_tensor[o2] = new RealMatrixBlock(num_block, left_index, right_index);

        delete[] left_index;
        delete[] right_index;
    }

    for(int o2=0;o2<right_bond;++o2)
    {
        for(int i=0;i<num_right_block;++i) delete[] zero_flag[o2][i];
        delete[] zero_flag[o2];
    }
    delete[] zero_flag;
    
    // compute left_contraction_tensor
    for(int o1=0;o1<left_bond;++o1)
    for(int i=0;i<left_contraction_tensor_[o1]->get_num_block();++i)
    {
        kdx = left_contracted_tensor_[o1]->left_index_[i];
        ldx = left_contracted_tensor_[o1]->right_index_[i];
        tmp_tensor[0] = left_contracted_tensor_[o1]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int k=0;k<num_same_index[kdx];++k)
        {
            idx = tensor_lattice->ket_tensor_->right_index_[position_same_index[kdx][k]];
            p1 = tensor_lattice->ket_tensor_->physics_index_[kdx][idx];
            if(tensor_operator->LeftCheckZero(o1, p1) == false) continue; //zero
            // L*B (b_0, a_0)'*(b_0, b1) = (a_0, b_1)
            tmp_tensor[1] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[kdx][k]];
            bra_tensor = tmp_tensor[1]->TransposeMatrix();
            first_tensor = bra_tensor->MultiplyToMatrix(tmp_tensor[0]);
            delete bra_tensor;
            for(int l=0;l<num_same_index[ldx];++l)
            {
                jdx = tensor_lattice->ket_tensor_->right_index_[position_same_index[ldx][l]];
                p2 = tensor_lattice->ket_tensor_->physics_index_[ldx][jdx];
                // L*B*K (a_0, b_1)*(b_1, a_1) = (a_0, a_1)
                tmp_tensor[2] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[ldx][l]];
                second_tensor = first_tensor->MultiplyToMatrix(tmp_tensor[2]);
                for(int o2=0;o2<right_bond;++o2)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_element(p1, p2);
                    if(element != 0.0) 
                    {
                        // L*B*K*W
                        second_tensor->MultiplyToScalar(element);
                        position = tensor_contraction[o2]->FindMatrixBlock(idx, jdx);
                        // matrix_block == nullptr?
                        tensor_contraction->left_contraction_tensor_[o2]->AddToMatrixBlock(position, second_tensor);
                    }
                }
                delete second_tensor;
            }
            delete first_tensor;
        }
    }

    delete[] num_same_index;
    for(int i=0;i<num_left_block;++i)
        delete[] position_same_index[i];
    delete position_same_index;
}



/*
//R^{o2, l, k}_{b_1, b_0}          R^{o1, j, i}_{a_1, a_0}
//                                          
//              ______                           ___ --- ______
//idx:k, dim:b_0      |            idx:i, dim:a_0   - B -      |
//                    |                              ---       |
//                    |                               |p1      |
//              ______|     _____\               ___-----______|
//       idx:o_2      |          /        idx:o_1   - W -      |
//                    |                             -----      |      
//                    |                               |p2      |
//              ______|                          ___ --- ______|
//idx:l, dim:b_1                   idx:j, dim:a_1   - K -
//                                                   ---  
*/                                       
void RealTensorContraction::RightComputeTensorContraction(RealTensorLattice* tensor_lattice, 
        RealTensorOperator* tensor_operator, RealTensorContraction* tensor_contraction)
{
    RealMatrix *first_tensor, *second_tensor, *tmp_tensor[3], *bra_tensor;
    bool ***zero_flag;
    int *num_same_index, **position_same_index, num_block, *left_index, *right_index;
    int left_bond, right_bond, num_left_block, num_right_block;
    int p1, p2, idx, jdx, kdx, ldx, index;
    double element;

    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();
    
    zero_flag = new int** [left_bond];
    for(int o=0;o<left_bond;++o)
    {
        zero_flag[o] = new int* [num_left_block];
        for(int j=0;j<num_left_block;++j)
        {
            zero_flag[o][i] = new int[num_left_block];
            for(int i=0;i<num_left_block;++i)
                zero_flag[o][j][i] = false;
        }
    }

    // for each right index(block), find all possible left indiced(blocks)
    num_same_index = new int[num_right_block];
    position_same_index = new int* [num_right_block];
    for(int i=0;i<num_left_block;++i)
        tensor_lattice->ket_tensor_->FindMatrixBlock(num_same_index[i], position_same_index[i], 1, i);

    // if zero_flag[o1][j][i] is true, 
    // then right_contraction_tensor[o1]->find_matrix_block(j, i) is nonezero
    for(int o2=0;o2<right_bond;++o2)
    for(int i=0;i<right_contraction_tensor_[o2]->get_num_block();++i)
    {
        ldx = right_contracted_tensor_[o2]->left_index_[i];
        kdx = right_contracted_tensor_[o2]->right_index_[i];
        tmp_tensor[0] = right_contracted_tensor_[o2]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int l=0;l<num_same_index[ldx];++l)
        {
            jdx = tensor_lattice->ket_tensor_->left_index_[position_same_index[ldx][l]];
            p2 = tensor_lattice->ket_tensor_->physics_index_[jdx][ldx];
            if(tensor_operator->RightCheckZero(o2, p2) == false) continue; //zero
            for(int k=0;k<num_same_index[kdx];++k)
            {
                idx = tensor_lattice->ket_tensor_->left_index_[position_same_index[kdx][k]];
                p1 = tensor_lattice->ket_tensor_->physics_index_[idx][kdx];
                for(int o1=0;o1<left_bond;++o1)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_element(p1, p2);
                    if(element != 0.0) zero_flag[o1][jdx][idx] = true;
                }
            }
        }
    }

    // generate right_contraction_tensor accronding to zero_flag,
    // (num_block, left_index, right_index)
    for(int o1=0;o1<left_bond;++o1)
    {
        num_block = 0;
        for(int j=0;j<num_left_block;++j) for(int i=0;i<num_left_block;++i)
        {
            if(zero_flag[o1][j][i] == true)
            {
                num_block++;
            }
        }

        if(num_block == 0)
        {
            left_index = nullptr;
            right_inex = nullptr;
        }
        else
        {
            index = 0;
            left_index = new int[num_block];
            right_index = new int[num_block];
            for(int j=0;j<num_left_block;++j) for(int i=0;i<num_left_block;++i)
            {
                if(zero_flag[o1][j][i] == true) 
                {
                    left_index[index] = j;
                    right_index[index] = i;
                    index++;
                }
            }
        }

        // if left_bond > right_bond?
        delete tensor_contraction->right_contraction_tensor[o1];
        tensor_contraction->left_contraction_tensor[o1] = new RealMatrixBlock(num_block, left_index, right_index);

        delete[] left_index;
        delete[] right_index;
    }

    for(int o1=0;o1<left_bond;++o1)
    {
        for(int j=0;j<num_left_block;++j) delete[] zero_flag[o1][j];
        delete[] zero_flag[o1];
    }
    delete[] zero_flag;
    
    // compute right_contraction_tensor
    for(int o2=0;o2<right_bond;++o2)
    for(int i=0;i<right_contraction_tensor_[o2]->get_num_block();++i)
    {
        ldx = right_contracted_tensor_[o2]->left_index_[i];
        kdx = right_contracted_tensor_[o2]->right_index_[i];
        tmp_tensor[0] = right_contracted_tensor_[o2]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int l=0;l<num_same_index[ldx];++l)
        {
            jdx = tensor_lattice->ket_tensor_->left_index_[position_same_index[ldx][l]];
            p2 = tensor_lattice->ket_tensor_->physics_index_[jdx][ldx];
            if(tensor_operator->RightCheckZero(o2, p2) == false) continue; //zero
            // K*R (a_1, b_1)*(b_1, b_0) = (a_1, b_0)
            tmp_tensor[1] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[ldx][l]];
            first_tensor = tmp_tensor[1]->MultiplyToMatrix(tmp_tensor[0]);
            for(int k=0;k<num_same_index[kdx];++k)
            {
                idx = tensor_lattice->ket_tensor_->left_index_[position_same_index[kdx][k]];
                p1 = tensor_lattice->ket_tensor_->physics_index_[idx][kdx];
                // K*R*B (a_1, b_0)*(a_0, b_0)' = (a_1, a_0)
                tmp_tensor[2] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[kdx][k]];
                bra_tensor = tmp_tensor[2]->TransposeMatrix();
                second_tensor = first_tensor->MultiplyToMatrix(bra_tensor);
                delete bra_tensor;
                for(int o1=0;o1<left_bond;++o1)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_element(p1, p2);
                    if(element != 0.0) 
                    {
                        // K*R*B*W
                        second_tensor->MultiplyToScalar(element);
                        position = tensor_contraction[o1]->FindMatrixBlock(jdx, idx);
                        // if matrix is nullptr, replace it
                        tensor_contraction->right_contraction_tensor_[o1]->AddToMatrixBlock(position, second_tensor);
                    }
                }
                delete second_tensor;
            }
            delete first_tensor;
        }
    }

    delete[] num_same_index;
    for(int i=0;i<num_right_block;++i)
        delete[] position_same_index[i];
    delete position_same_index;
}

/*
// effictive hamiltonian:
// H^{p1,p2,i,j,k,l}_{a_0,b_0,a_1,b_1}
//                ______     ______
//               | i,a_0     j,b_0 |
//               |                 |
//               |        |p1      |
//               |______-----______|
//               |  o1  - W -  o2  |
//               |      -----      |      
//               |        |p2      |
//               |______     ______|
//                 k,a_1     l,b_1          
//                        
*/
void RealTensorContraction::ComputeEffectHamilton(TensorLattice* tensor_lattice, 
        TensorOperator* tensor_operator, double* hamilton)
{
    RealMatrix *inv_tensor, *result_tensor, *tmp_tensor[2];
    double element, *matrix_element;
    int vector_dim, position[2], result_dim[2], part_dim[2];
    int p1, p2, idx, jdx, kdx, ldx;

    vector_dim = tensor_lattice->ComputeKetTensorDim();

    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();
    for(int o1=0;o1<left_bond;++o1)
    for(int i=0;i<left_contraction_tensor_[o1]->get_num_block();++i)
    {
        idx = left_contraction_tensor_[o1]->left_index_[i];
        kdx = left_contraction_tensor_[o1]->right_index_[i];
        tmp_tensor[0] = left_contraction_tensor_[o1]->matrix_block_[i];

        if(tmp_tensor[0]->get_row() == 0) continue;

        for(int o2=0;o2<right_bond;++o2)
        for(int j=0;j<right_contraction_tensor_[o2]->get_num_block();++j)
        {
            ldx = right_contraction_tensor_[o2]->left_index_[j];
            jdx = right_contraction_tensor_[o2]->right_index_[j];
            tmp_tensor[1] = right_contraction_tensor_[o2]->matrix_block_[j];

            if(tmp_tensor[1]->get_row() == 0) continue;
            
            // if exit (p1,i,j) and (p2,k,l)
            position[0] = tensor_lattice->FindMatrixBlock(idx, jdx);
            position[1] = tensor_lattice->FindMatrixBlock(kdx, ldx);

            if(position[0]!=-1 && position[1]!=-1)
            {
                p1 = tensor_lattice->physics_index_[idx][jdx];
                p2 = tensor_lattice->physics_index_[kdx][ldx];
                element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                if(element != 0)
                {
                    inv_tensor = tmp_tensor[1]->TransposeMatrix();
                    result_tensor = tmp_tensor[0]->MultiplyToMatrix(inv_tensor);
                    delete inv_tensor;

                    if(result_tensor != nullptr)
                    {
                        result_tensor->MultiplyToScale(element);

                        matrix_element = result_tensor->get_matrix_element();
                        
                        result_dim[0] = result_tensor->get_row();
                        result_dim[1] = result_tensor->get_column();

                        part_dim[0] = tensor_lattice->ComputePartKetTensorDim(position[0]);
                        part_dim[1] = tensor_lattice->ComputePartKetTensorDim(position[1]);

                        for(int a=0;a<result_dim[0];++a) for(int b=0;b<result_dim[1];++b)
                        {
                            // hamilton is zero?
                            hamilton[part_dim[1]+b + (part_dim[0]+a)*vector_dim] += matrix_element[b+a*result_dim[1]];
                        }
                    }
                    delete result_tensor;
                }
            }
        }
    }
}


/*
// effictive hamiltonian multiply to state 
// B^{p1,i,j}_{a_0,b_0} = L*W*K*R
//                ______     ______
//               | i,a_0     j,b_0 |
//               |                 |
//               |        |p1      |
//               |______-----______|
//               |  o1  - W -  o2  |
//               |      -----      |      
//               |        |p2      |
//               |______ --- ______|
//                 k,a_1- K -l,b_1          
//                       --- 
*/
void RealTensorContraction::MultiplyEffectHamilton(RealTensorLattice* tensor_lattice, 
        RealTensorOperator* tensor_operator, double* state)
{
    RealMatrix *first_matrix, *second_matrix, tmp_tensor[3];
    double *matrix_element, element;
    int *num_same_index, **position_same_index, part_dim, position, total_element_num, p1, p2, idx, jdx, ldx, kdx;


    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();

    // for each right index(block), find all possible left indiced(blocks)
    num_same_index = new int[num_right_block];
    position_same_index = new int* [num_right_block];
    for(int i=0;i<num_left_block;++i)
        tensor_lattice->ket_tensor_->FindMatrixBlock(num_same_index[i], position_same_index[i], 1, i);

    // compute right_contraction_tensor
    for(int o2=0;o2<right_bond;++o2)
    for(int i=0;i<right_contraction_tensor_[o2]->get_num_block();++i)
    {
        ldx = right_contracted_tensor_[o2]->left_index_[i];
        jdx = right_contracted_tensor_[o2]->right_index_[i];
        tmp_tensor[0] = right_contracted_tensor_[o2]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int l=0;l<num_same_index[ldx];++l)
        {
            kdx = tensor_lattice->ket_tensor_->left_index_[position_same_index[ldx][l]];
            p2 = tensor_lattice->ket_tensor_->physics_index_[kdx][ldx];
            if(tensor_operator->RightCheckZero(o2, p2) == false) continue; //zero
            // K*R (a_1, b_1)*(b_1, b_0) = (a_1, b_0)
            tmp_tensor[1] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[ldx][l]];
            first_tensor = tmp_tensor[1]->MultiplyToMatrix(tmp_tensor[0]);
            for(int j=0;j<num_same_index[jdx];++j)
            {
                idx = tensor_lattice->ket_tensor_->left_index_[position_same_index[jdx][k]];
                p1 = tensor_lattice->ket_tensor_->physics_index_[idx][jdx];
                part_dim = tensor_lattice->ComputePartKetTensorDim(position_same_index[jdx][k])
                for(int o1=0;o1<left_bond;++o1)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_element(p1, p2);
                    if(element != 0.0) 
                    {
                        position = tensor_contraction[o1]->FindMatrixBlock(idx, kdx);
                        if(position != -1)
                        {
                            tmp_tensor[2] = left_contraction_tensor[o1]->matrix_block_[position];
                            second_tensor = tmp_tensor[2]->MultiplyToMatrix(first_tensor);
                            second_tensor->MultiplyToScale(element);
                            matrix_element = second_tensor->get_matrix_element();
                            total_element_num = second_tensor->get_total_element_num();
                            for(int a=0;a<total_element_num;++a) state[part_dim+a] = matrix_element[a];
                            delete second_tensor;
                        }
                    }
                }
            }
            delete first_tensor;
        }
    }

    delete[] num_same_index;
    for(int i=0;i<num_right_block;++i)
        delete[] position_same_index[i];
    delete position_same_index;
}
