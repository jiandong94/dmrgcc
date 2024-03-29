#include "dmrg/complexdmrg/complex_tensor_contraction.h"

ComplexTensorContraction::ComplexTensorContraction(int left_bond, int right_bond)
{
    left_bond_ = left_bond;
    right_bond_ = right_bond;

    left_contraction_tensor_ = new ComplexMatrixBlock* [left_bond_];
    for(int i=0;i<left_bond_;++i)
        left_contraction_tensor_[i] = nullptr;

    right_contraction_tensor_ = new ComplexMatrixBlock* [right_bond_];
    for(int i=0;i<right_bond_;++i)
        right_contraction_tensor_[i] = nullptr;
}

ComplexTensorContraction::~ComplexTensorContraction()
{
    for(int i=0;i<left_bond_;++i)
        delete left_contraction_tensor_[i];
    delete[] left_contraction_tensor_;
    
    for(int i=0;i<right_bond_;++i)
        delete right_contraction_tensor_[i];
    delete[] right_contraction_tensor_;
}

int ComplexTensorContraction::get_left_bond()
{
    return left_bond_;
}

int ComplexTensorContraction::get_right_bond()
{
    return right_bond_;
}

ComplexMatrixBlock** ComplexTensorContraction::get_left_contraction_tensor()
{
    return left_contraction_tensor_;
}

ComplexMatrixBlock** ComplexTensorContraction::get_right_contraction_tensor()
{
    return right_contraction_tensor_;
}

void ComplexTensorContraction::PrintTensorContraction()
{
    cout << "left bond = " << left_bond_ << ", ";
    cout << "right bond = " << right_bond_ << endl;
    cout << " ----------------------------" << endl;
    cout << "|left contraction tensor L : | " << endl;
    cout << " ----------------------------" << endl;
    for(int i=0;i<left_bond_;++i)
    {
        cout << "left bond index = " << i << endl;
        if(left_contraction_tensor_[i] != nullptr)
            left_contraction_tensor_[i]->PrintMatrixBlock();
    }
    
    cout << " ----------------------------" << endl; 
    cout << "|right contraction tensor R :| " << endl;
    cout << " ----------------------------" << endl;
    for(int i=0;i<right_bond_;++i)
    {
        cout << "right bond index = " << i << endl;
        if(right_contraction_tensor_[i] != nullptr)
            right_contraction_tensor_[i]->PrintMatrixBlock();
    }
}

void ComplexTensorContraction::WriteTensorContraction(int leigh, char* tensor_contraction_name)
{
    ofstream tensor_contraction_file;
    tensor_contraction_file.open(tensor_contraction_name, ios::binary|ios::out);

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

void ComplexTensorContraction::ReadTensorContraction(int leigh, char* tensor_contraction_name)
{
    ifstream tensor_contraction_file;
    tensor_contraction_file.open(tensor_contraction_name, ios::binary|ios::in);

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

void ComplexTensorContraction::DefineTensorContraction(int leigh)
{
    ComplexMatrix *tmp_tensor;
    
    tmp_tensor = new ComplexMatrix(1,1);
    tmp_tensor->set_matrix_element(0,0,1);

    if(leigh == 0)
    {
        for(int i=0;i<left_bond_;++i)
        {
            if(left_contraction_tensor_[i] != nullptr) delete left_contraction_tensor_[i];
            left_contraction_tensor_[i] = new ComplexMatrixBlock(1);
            left_contraction_tensor_[i]->set_matrix_block(0, tmp_tensor);
        }
    }
    else if(leigh == 1)
    {
        for(int i=0;i<right_bond_;++i)
        {
            if(right_contraction_tensor_[i] != nullptr) delete right_contraction_tensor_[i];
            right_contraction_tensor_[i] = new ComplexMatrixBlock(1);
            right_contraction_tensor_[i]->set_matrix_block(0, tmp_tensor);
        }
    }
    delete tmp_tensor;
}

void ComplexTensorContraction::ResetTensorContraction(int leigh)
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
//              E^{p1,i,o2+a_1}_{a_0,a_1}   
//                                          
//                ______ 
//               |i,a_0
//               |      
//               |        |p1
//               |______-----___      __  
//               | o_1  - W -   o_2     |
//               |      -----           |_ o_2+j 
//               |        |p2           |
//               |______ --- ___      __|
//                k,b_1 - K -   j,a_1
//                       ---  
*/                                       
void ComplexTensorContraction::LeftExpanTensorContraction(ComplexTensorLattice* tensor_lattice, ComplexTensorOperator* tensor_operator, 
        ComplexTensorLattice* expan_tensor_lattice, int** mapping_table, double noise_factor)
{
    ComplexMatrix *expan_tensor, *first_tensor, *second_tensor, *tmp_tensor[3];
    bool *exist_flag;
    Complex element;
    int left_bond, right_bond, num_left_block, num_right_block;
    int *num_same_index, **position_same_index, **expan_table, min_expan, max_expan;
    int num_block, *left_index, *right_index, **physics_index;
    int index, first_dim[2], expan_dim[2], num_leigh_block, *leigh_block;
    int position, p1, p2, idx, jdx, kdx, ldx;

    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();
    
    exist_flag = new bool[expan_tensor_lattice->get_num_right_block()];
    for(int i=0;i<expan_tensor_lattice->get_num_right_block();++i)
        exist_flag[i] = false;
    
    // for each left index(block), find all possible right indiced(blocks)
    num_same_index = new int[num_left_block];
    position_same_index = new int* [num_left_block];
    for(int i=0;i<num_left_block;++i)
        tensor_lattice->ket_tensor_->FindMatrixBlock(num_same_index[i], position_same_index[i], 0, i);

    // expan_table is like mapping_table, which 
    expan_table = new int* [right_bond];
    for(int o2=0;o2<right_bond;++o2)
    {
        expan_table[o2] = new int[num_right_block];
        for(int r=0;r<num_right_block;++r)
            expan_table[o2][r] = -1;
    }

    min_expan = 1;
    max_expan = 25;
    //--------------------------------------------------------------------------------------
    // exist_flag
    for(int o1=0;o1<left_bond;++o1)
    for(int i=0;i<left_contraction_tensor_[o1]->get_num_block();++i)
    {
        idx = left_contraction_tensor_[o1]->left_index_[i];
        kdx = left_contraction_tensor_[o1]->right_index_[i];
        tmp_tensor[0] = left_contraction_tensor_[o1]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int k=0;k<num_same_index[kdx];++k)
        {
            jdx = tensor_lattice->ket_tensor_->right_index_[position_same_index[kdx][k]];
            p2 = tensor_lattice->physics_index_[kdx][jdx];
            if(tensor_operator->MixedCheckZero(o1, p2) == false) continue; //zero
            for(int o2=0;o2<right_bond;++o2)
            {
                ldx = mapping_table[o2][jdx];
                position = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(idx, ldx);

                if(ldx!=-1 && position!=-1)
                {
                    p1 = expan_tensor_lattice->physics_index_[idx][ldx];
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0) exist_flag[ldx] = true;
                }
            }
        }
    }

    // num_leigh_block, leigh_block, physics_index
    num_leigh_block = 0;
    for(int i=0;i<expan_tensor_lattice->get_num_right_block();++i)
        if(exist_flag[i] != false)
            num_leigh_block++;
    leigh_block = new int[num_leigh_block];
    physics_index = new int* [expan_tensor_lattice->num_left_block_];
    for(int i=0;i<expan_tensor_lattice->num_left_block_;++i)
    {
        physics_index[i] = new int[num_leigh_block];
        for(int j=0;j<num_leigh_block;++j)
            physics_index[i][j] = -1;
    }

    index = 0;
    for(int i=0;i<expan_tensor_lattice->get_num_right_block();++i)
    {
        if(exist_flag[i] == true)
        {
            leigh_block[index] = expan_tensor_lattice->right_block_[i];
            for(int j=0;j<expan_tensor_lattice->get_num_left_block();++j)
            {
                physics_index[j][index] = expan_tensor_lattice->physics_index_[j][i];
            }
            index++;
        }
    }

    // expan_table
    for(int o2=0;o2<right_bond;++o2) for(int r=0;r<num_right_block;++r)
    {
        expan_table[o2][r] = -1;
        ldx = mapping_table[o2][r];
        if(ldx != -1)
        {
            for(int i=0;i<num_leigh_block;++i)
            {
                if(expan_tensor_lattice->right_block_[ldx] == leigh_block[i])
                {
                    expan_table[o2][r] = i;
                    break;
                }
            }
        }

    }
    
    // re-construct left_dim, right_block, right_dim, physics_index
    // for expan_lattice
    delete[] expan_tensor_lattice->right_block_;
    delete[] expan_tensor_lattice->right_dim_;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i) 
        delete[] expan_tensor_lattice->physics_index_[i];
    delete[] expan_tensor_lattice->physics_index_;

    expan_tensor_lattice->num_right_block_ = num_leigh_block;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
        expan_tensor_lattice->left_dim_[i] = tensor_lattice->left_dim_[i];
    expan_tensor_lattice->right_block_ = new int[expan_tensor_lattice->num_right_block_];
    expan_tensor_lattice->right_dim_ = new int[expan_tensor_lattice->num_right_block_];
    for(int i=0;i<expan_tensor_lattice->num_right_block_;++i)
    {
        expan_tensor_lattice->right_block_[i] = leigh_block[i];
        expan_tensor_lattice->right_dim_[i] = 0;
    }
    
    expan_tensor_lattice->physics_index_ = new int* [expan_tensor_lattice->num_left_block_];
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
    {
        expan_tensor_lattice->physics_index_[i] = new int[expan_tensor_lattice->num_right_block_];
        for(int j=0;j<expan_tensor_lattice->get_num_right_block();++j)
        {
            expan_tensor_lattice->physics_index_[i][j] = physics_index[i][j];
        }
    }
    
    // re-construct ket_tensor for expan_lattice
    delete expan_tensor_lattice->ket_tensor_;
    num_block = 0;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
    for(int j=0;j<expan_tensor_lattice->get_num_right_block();++j)
    {
        if(expan_tensor_lattice->physics_index_[i][j] != -1)
            num_block++;
    }

    left_index = nullptr;
    right_index = nullptr;
    if(num_block != 0)
    {
        left_index = new int[num_block];
        right_index = new int[num_block];
        index = 0;
        for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
        for(int j=0;j<expan_tensor_lattice->get_num_right_block();++j)
        {
            if(expan_tensor_lattice->physics_index_[i][j] != -1)
            {
                left_index[index] = i;
                right_index[index] = j;
                index++;
            }
        }
        
    }
    expan_tensor_lattice->ket_tensor_ = new ComplexMatrixBlock(num_block, left_index, right_index);
    delete[] exist_flag;
    delete[] leigh_block;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
        delete[] physics_index[i];
    delete[] physics_index;
    delete[] left_index;
    delete[] right_index;

    //--------------------------------------------------------------------------------------
    //
    for(int o1=0;o1<left_bond;++o1)
    for(int i=0;i<left_contraction_tensor_[o1]->get_num_block();++i)
    {
        idx = left_contraction_tensor_[o1]->left_index_[i];
        kdx = left_contraction_tensor_[o1]->right_index_[i];
        tmp_tensor[0] = left_contraction_tensor_[o1]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int k=0;k<num_same_index[kdx];++k)
        {
            jdx = tensor_lattice->ket_tensor_->right_index_[position_same_index[kdx][k]];
            p2 = tensor_lattice->physics_index_[kdx][jdx];
            if(tensor_operator->MixedCheckZero(o1, p2) == false) continue; //zero

            tmp_tensor[1] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[kdx][k]];
            
            // L*K (a_0*b_1)*(b_1*a_1) = (a_0*a_1)
            first_tensor = tmp_tensor[0]->MultiplyToMatrix(tmp_tensor[1]);
            
            first_dim[0] = first_tensor->get_row();
            first_dim[1] = first_tensor->get_column();

            expan_dim[0] = first_dim[1];
            expan_dim[1] = (int)(sqrt(first_dim[1]));
            expan_dim[1] *= tensor_lattice->physics_dim_;
            expan_dim[1] *= tensor_lattice->physics_dim_;
            if(expan_dim[1] < min_expan) expan_dim[1] = min_expan;
            if(expan_dim[1] > max_expan) expan_dim[1] = max_expan;

            tmp_tensor[2] = new ComplexMatrix(expan_dim[0], expan_dim[1]);
            //cout << "Liden: " << min(expan_dim[0], expan_dim[1]) << endl;
            //cout << "total element: " << tmp_tensor[2]->get_total_element_num() << endl;
            
            for(int o2=0;o2<right_bond;++o2)
            {
                ldx = expan_table[o2][jdx];
                position = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(idx, ldx);
                
                if(ldx!=-1 && position!=-1)
                {
                    p1 = expan_tensor_lattice->physics_index_[idx][ldx];
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0)
                    {
                        tmp_tensor[2]->RandomMatrix();
                        for(int t=0;t<min(expan_dim[0], expan_dim[1]);++t)
                            tmp_tensor[2]->set_matrix_element(t, t, 0.1*element);
                        // L*K*W*random_matrix
                        second_tensor = first_tensor->MultiplyToMatrix(tmp_tensor[2]);
                        expan_tensor = expan_tensor_lattice->ket_tensor_->matrix_block_[position];
                        if(expan_tensor->get_row() == 0)
                            expan_tensor->AddToMatrix(second_tensor);
                        else
                            expan_tensor->ExpanMatrix(1, second_tensor);
                        delete second_tensor;
                    }

                }
            }
            delete first_tensor;
            delete tmp_tensor[2];
        }
    }
    delete[] num_same_index;
    for(int i=0;i<num_left_block;++i)
        delete[] position_same_index[i];
    delete position_same_index;

    for(int i=0;i<right_bond;++i)
        delete[] expan_table[i];
    delete[] expan_table;
    //---------------------------------------------------------------------------------------
    // re-compute right_dim
    for(int r=0;r<expan_tensor_lattice->get_num_right_block();++r)
    {
        expan_tensor_lattice->right_dim_[r] = 0;
        for(int l=0;l<expan_tensor_lattice->get_num_left_block();++l)
        {
            position = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(l, r);
            if(position != -1)
            {
                expan_tensor = expan_tensor_lattice->ket_tensor_->matrix_block_[position];
                expan_dim[1] = expan_tensor->get_column();
                
                if(expan_dim[1] > expan_tensor_lattice->right_dim_[r]) 
                    expan_tensor_lattice->right_dim_[r] = expan_dim[1];
            }
        }
        
        for(int l=0;l<expan_tensor_lattice->get_num_left_block();++l)
        {
            position = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(l, r);
            if(position != -1)
            {
                expan_tensor = expan_tensor_lattice->ket_tensor_->matrix_block_[position];
                if(expan_tensor->get_row() == 0)
                {
                    cout << "row of expan_tensor is zero in LeftExpanTensorContraction" << endl;
                    first_tensor = new ComplexMatrix(expan_tensor_lattice->left_dim_[l], expan_tensor_lattice->right_dim_[r]);
                    expan_tensor->AddToMatrix(first_tensor);
                    delete first_tensor;
                }
                else
                {
                    expan_dim[1] = expan_tensor_lattice->right_dim_[r] - expan_tensor->get_column();
                    if(expan_dim[1] > 0)
                    {
                        first_tensor = new ComplexMatrix(expan_tensor->get_row(), expan_dim[1]);
                        expan_tensor->ExpanMatrix(1, first_tensor);
                        delete first_tensor;
                    }
                }
            }
        }
    } 

    expan_tensor_lattice->ket_tensor_->MultiplyToScalar(noise_factor);
}


/*
//                    F^{p1,o1+j,i}_{a_1,a_0}   
//                                          
//                                    ______
//                                     i,a_0|
//                                          |
//                                 |p1      |
//                     _      ___-----______|       
//                    |    o_1   - W -  o_2 |   
//              o_1+j_|          -----      |     
//                    |            |p2      |     
//                    |_      ___ --- ______|      
//                       j,a_1   - K - k,b_1 
//                                ---  
*/                                       
void ComplexTensorContraction::RightExpanTensorContraction(ComplexTensorLattice* tensor_lattice, ComplexTensorOperator* tensor_operator, 
        ComplexTensorLattice* expan_tensor_lattice, int** mapping_table, double noise_factor)
{
    ComplexMatrix *expan_tensor, *first_tensor, *second_tensor, *tmp_tensor[3];
    bool *exist_flag;
    Complex element;
    int left_bond, right_bond, num_left_block, num_right_block;
    int *num_same_index, **position_same_index, **expan_table, min_expan, max_expan;
    int num_block, *left_index, *right_index, **physics_index;
    int index, first_dim[2], expan_dim[2], num_leigh_block, *leigh_block;
    int position, p1, p2, idx, jdx, kdx, ldx;

    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();
    
    exist_flag = new bool[expan_tensor_lattice->get_num_left_block()];
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
        exist_flag[i] = false;
    
    // for each right index(block), find all possible left indiced(blocks)
    num_same_index = new int[num_right_block];
    position_same_index = new int* [num_right_block];
    for(int i=0;i<num_right_block;++i)
        tensor_lattice->ket_tensor_->FindMatrixBlock(num_same_index[i], position_same_index[i], 1, i);

    // expan_table is like mapping_table, which 
    expan_table = new int* [left_bond];
    for(int o1=0;o1<left_bond;++o1)
    {
        expan_table[o1] = new int[num_left_block];
        for(int l=0;l<num_left_block;++l)
            expan_table[o1][l] = -1;
    }

    min_expan = 1;
    max_expan = 25;
    //--------------------------------------------------------------------------------------
    // exist_flag
    for(int o2=0;o2<right_bond;++o2)
    for(int i=0;i<right_contraction_tensor_[o2]->get_num_block();++i)
    {
        kdx = right_contraction_tensor_[o2]->left_index_[i];
        idx = right_contraction_tensor_[o2]->right_index_[i];
        tmp_tensor[0] = right_contraction_tensor_[o2]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int k=0;k<num_same_index[kdx];++k)
        {
            jdx = tensor_lattice->ket_tensor_->left_index_[position_same_index[kdx][k]];
            p2 = tensor_lattice->physics_index_[jdx][kdx];
            if(tensor_operator->RightCheckZero(o2, p2) == false) continue; //zero
            for(int o1=0;o1<left_bond;++o1)
            {
                ldx = mapping_table[o1][jdx];
                position = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(ldx, idx);

                if(ldx!=-1 && position!=-1)
                {
                    p1 = expan_tensor_lattice->physics_index_[ldx][idx];
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0) exist_flag[ldx] = true;
                }
            }
        }
    }

    // num_leigh_block, leigh_block, physics_index
    num_leigh_block = 0;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
        if(exist_flag[i] != false)
            num_leigh_block++;
    leigh_block = new int[num_leigh_block];
    physics_index = new int* [num_leigh_block];
    for(int i=0;i<num_leigh_block;++i)
    {
        physics_index[i] = new int[expan_tensor_lattice->num_right_block_];
        for(int j=0;j<expan_tensor_lattice->num_right_block_;++j)
            physics_index[i][j] = -1;
    }

    index = 0;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
    {
        if(exist_flag[i] == true)
        {
            leigh_block[index] = expan_tensor_lattice->left_block_[i];
            for(int j=0;j<expan_tensor_lattice->get_num_right_block();++j)
            {
                physics_index[index][j] = expan_tensor_lattice->physics_index_[i][j];
            }
            index++;
        }
    }

    // expan_table
    for(int o1=0;o1<left_bond;++o1) for(int l=0;l<num_left_block;++l)
    {
        expan_table[o1][l] = -1;
        ldx = mapping_table[o1][l];
        if(ldx != -1)
        {
            for(int i=0;i<num_leigh_block;++i)
            {
                if(expan_tensor_lattice->left_block_[ldx] == leigh_block[i])
                {
                    expan_table[o1][l] = i;
                    break;
                }
            }
        }

    }
    
    // re-construct left_dim, right_block, right_dim, physics_index
    // for expan_lattice
    delete[] expan_tensor_lattice->left_block_;
    delete[] expan_tensor_lattice->left_dim_;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i) 
        delete[] expan_tensor_lattice->physics_index_[i];
    delete[] expan_tensor_lattice->physics_index_;

    expan_tensor_lattice->num_left_block_ = num_leigh_block;
    for(int i=0;i<expan_tensor_lattice->get_num_right_block();++i)
        expan_tensor_lattice->right_dim_[i] = tensor_lattice->right_dim_[i];
    expan_tensor_lattice->left_block_ = new int[expan_tensor_lattice->num_left_block_];
    expan_tensor_lattice->left_dim_ = new int[expan_tensor_lattice->num_left_block_];
    for(int i=0;i<expan_tensor_lattice->num_left_block_;++i)
    {
        expan_tensor_lattice->left_block_[i] = leigh_block[i];
        expan_tensor_lattice->left_dim_[i] = 0;
    }
    
    expan_tensor_lattice->physics_index_ = new int* [expan_tensor_lattice->num_left_block_];
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
    {
        expan_tensor_lattice->physics_index_[i] = new int[expan_tensor_lattice->num_right_block_];
        for(int j=0;j<expan_tensor_lattice->get_num_right_block();++j)
        {
            expan_tensor_lattice->physics_index_[i][j] = physics_index[i][j];
        }
    }

    // re-construct ket_tensor for expan_lattice
    delete expan_tensor_lattice->ket_tensor_;
    num_block = 0;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
    for(int j=0;j<expan_tensor_lattice->get_num_right_block();++j)
    {
        if(expan_tensor_lattice->physics_index_[i][j] != -1)
            num_block++;
    }

    left_index = nullptr;
    right_index = nullptr;
    if(num_block != 0)
    {
        left_index = new int[num_block];
        right_index = new int[num_block];
        index = 0;
        for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
        for(int j=0;j<expan_tensor_lattice->get_num_right_block();++j)
        {
            if(expan_tensor_lattice->physics_index_[i][j] != -1)
            {
                left_index[index] = i;
                right_index[index] = j;
                index++;
            }
        }
        
    }
    expan_tensor_lattice->ket_tensor_ = new ComplexMatrixBlock(num_block, left_index, right_index);

    delete[] exist_flag;
    delete[] leigh_block;
    for(int i=0;i<expan_tensor_lattice->get_num_left_block();++i)
        delete[] physics_index[i];
    delete[] physics_index;
    delete[] left_index;
    delete[] right_index;

    //--------------------------------------------------------------------------------------
    //
    for(int o2=0;o2<right_bond;++o2)
    for(int i=0;i<right_contraction_tensor_[o2]->get_num_block();++i)
    {
        kdx = right_contraction_tensor_[o2]->left_index_[i];
        idx = right_contraction_tensor_[o2]->right_index_[i];
        tmp_tensor[0] = right_contraction_tensor_[o2]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int k=0;k<num_same_index[kdx];++k)
        {
            jdx = tensor_lattice->ket_tensor_->left_index_[position_same_index[kdx][k]];
            p2 = tensor_lattice->physics_index_[jdx][kdx];
            if(tensor_operator->RightCheckZero(o2, p2) == false) continue; //zero

            tmp_tensor[1] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[kdx][k]];
            
            // K*R (a_1, b_1)*(b_1, a_0) = (a_1, a_0)
            first_tensor = tmp_tensor[1]->MultiplyToMatrix(tmp_tensor[0]);
            first_dim[0] = first_tensor->get_row();
            first_dim[1] = first_tensor->get_column();

            expan_dim[1] = first_dim[0];
            expan_dim[0] = (int)(sqrt(first_dim[0]));
            expan_dim[0] *= tensor_lattice->physics_dim_;
            expan_dim[0] *= tensor_lattice->physics_dim_;
            if(expan_dim[0] < min_expan) expan_dim[0] = min_expan;
            if(expan_dim[0] > max_expan) expan_dim[0] = max_expan;

            tmp_tensor[2] = new ComplexMatrix(expan_dim[0], expan_dim[1]);

            for(int o1=0;o1<left_bond;++o1)
            {
                ldx = expan_table[o1][jdx];
                position = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(ldx, idx);

                if(ldx!=-1 && position!=-1)
                {
                    p1 = expan_tensor_lattice->physics_index_[ldx][idx];
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0)
                    {
                        tmp_tensor[2]->RandomMatrix();
                        for(int t=0;t<min(expan_dim[0], expan_dim[1]);++t)
                            tmp_tensor[2]->set_matrix_element(t, t, 0.1*element);
                        // K*R*W*random_matrix
                        second_tensor = tmp_tensor[2]->MultiplyToMatrix(first_tensor);
                        expan_tensor = expan_tensor_lattice->ket_tensor_->matrix_block_[position];
                        if(expan_tensor->get_row() == 0)
                            expan_tensor->AddToMatrix(second_tensor);
                        else 
                            expan_tensor->ExpanMatrix(0, second_tensor);
                        delete second_tensor;
                    }
                }
            }
            delete first_tensor;
            delete tmp_tensor[2];
        }
    }
    delete[] num_same_index;
    for(int i=0;i<num_right_block;++i)
        delete[] position_same_index[i];
    delete position_same_index;

    for(int i=0;i<left_bond;++i)
        delete[] expan_table[i];
    delete[] expan_table;

    //---------------------------------------------------------------------------------------
    // re-compute left_dim
    for(int l=0;l<expan_tensor_lattice->get_num_left_block();++l)
    {
        expan_tensor_lattice->left_dim_[l] = 0;
        for(int r=0;r<expan_tensor_lattice->get_num_right_block();++r)
        {
            position = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(l, r);
            if(position != -1)
            {
                expan_tensor = expan_tensor_lattice->ket_tensor_->matrix_block_[position];
                expan_dim[0] = expan_tensor->get_row();
                if(expan_dim[0] > expan_tensor_lattice->left_dim_[l]) 
                    expan_tensor_lattice->left_dim_[l] = expan_dim[0];
            }
        }
        
        for(int r=0;r<expan_tensor_lattice->get_num_right_block();++r)
        {
            position = expan_tensor_lattice->ket_tensor_->FindMatrixBlock(l, r);
            if(position != -1)
            {
                expan_tensor = expan_tensor_lattice->ket_tensor_->matrix_block_[position];
                if(expan_tensor->get_row() == 0)
                {
                    cout << "row of expan_tensor is zero in RightExpanTensorContraction" << endl;
                    first_tensor = new ComplexMatrix(expan_tensor_lattice->left_dim_[l], expan_tensor_lattice->right_dim_[r]);
                    expan_tensor->AddToMatrix(first_tensor);
                    delete first_tensor;
                }
                else
                {
                    expan_dim[0] = expan_tensor_lattice->left_dim_[l] - expan_tensor->get_row();
                    if(expan_dim[0] > 0)
                    {
                        first_tensor = new ComplexMatrix(expan_dim[0], expan_tensor->get_column());
                        expan_tensor->ExpanMatrix(0, first_tensor);
                        delete first_tensor;
                    }
                }
            }
        }
    } 

    expan_tensor_lattice->ket_tensor_->MultiplyToScalar(noise_factor);
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
void ComplexTensorContraction::LeftComputeTensorContraction(ComplexTensorLattice* tensor_lattice, 
        ComplexTensorOperator* tensor_operator, ComplexTensorContraction* tensor_contraction)
{
    ComplexMatrix *first_tensor, *second_tensor, *tmp_tensor[3], *bra_tensor;
    bool ***zero_flag;
    int *num_same_index, **position_same_index, num_block, *left_index, *right_index;
    int left_bond, right_bond, num_left_block, num_right_block;
    int p1, p2, idx, jdx, kdx, ldx, index, position;
    Complex element;
    
    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();
    
    zero_flag = new bool** [right_bond];
    for(int o=0;o<right_bond;++o)
    {
        zero_flag[o] = new bool* [num_right_block];
        for(int i=0;i<num_right_block;++i)
        {
            zero_flag[o][i] = new bool[num_right_block];
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
        kdx = left_contraction_tensor_[o1]->left_index_[i];
        ldx = left_contraction_tensor_[o1]->right_index_[i];
        tmp_tensor[0] = left_contraction_tensor_[o1]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int k=0;k<num_same_index[kdx];++k)
        {
            idx = tensor_lattice->ket_tensor_->right_index_[position_same_index[kdx][k]];
            p1 = tensor_lattice->physics_index_[kdx][idx];
            if(tensor_operator->LeftCheckZero(o1, p1) == false) continue; //zero
            for(int l=0;l<num_same_index[ldx];++l)
            {
                jdx = tensor_lattice->ket_tensor_->right_index_[position_same_index[ldx][l]];
                p2 = tensor_lattice->physics_index_[ldx][jdx];
                for(int o2=0;o2<right_bond;++o2)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0) zero_flag[o2][idx][jdx] = true;
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
            right_index = nullptr;
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
        delete tensor_contraction->left_contraction_tensor_[o2];
        tensor_contraction->left_contraction_tensor_[o2] = new ComplexMatrixBlock(num_block, left_index, right_index);

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
        kdx = left_contraction_tensor_[o1]->left_index_[i];
        ldx = left_contraction_tensor_[o1]->right_index_[i];
        tmp_tensor[0] = left_contraction_tensor_[o1]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int k=0;k<num_same_index[kdx];++k)
        {
            idx = tensor_lattice->ket_tensor_->right_index_[position_same_index[kdx][k]];
            p1 = tensor_lattice->physics_index_[kdx][idx];
            if(tensor_operator->LeftCheckZero(o1, p1) == false) continue; //zero
            // L*B (b_0, a_0)'*(b_0, b1) = (a_0, b_1)
            tmp_tensor[1] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[kdx][k]];
            bra_tensor = tmp_tensor[1]->HermitianConjugateMatrix();
            first_tensor = bra_tensor->MultiplyToMatrix(tmp_tensor[0]);
            delete bra_tensor;
            for(int l=0;l<num_same_index[ldx];++l)
            {
                jdx = tensor_lattice->ket_tensor_->right_index_[position_same_index[ldx][l]];
                p2 = tensor_lattice->physics_index_[ldx][jdx];
                // L*B*K (a_0, b_1)*(b_1, a_1) = (a_0, a_1)
                tmp_tensor[2] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[ldx][l]];
                second_tensor = first_tensor->MultiplyToMatrix(tmp_tensor[2]);
                for(int o2=0;o2<right_bond;++o2)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0) 
                    {
                        // L*B*K*W
                        position = tensor_contraction->left_contraction_tensor_[o2]->FindMatrixBlock(idx, jdx);
                        // matrix_block == nullptr?
                        tensor_contraction->left_contraction_tensor_[o2]->
                                            AddToMatrixBlock(position, element, second_tensor);
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
void ComplexTensorContraction::RightComputeTensorContraction(ComplexTensorLattice* tensor_lattice, 
        ComplexTensorOperator* tensor_operator, ComplexTensorContraction* tensor_contraction)
{
    ComplexMatrix *first_tensor, *second_tensor, *tmp_tensor[3], *bra_tensor;
    bool ***zero_flag;
    int *num_same_index, **position_same_index, num_block, *left_index, *right_index;
    int left_bond, right_bond, num_left_block, num_right_block;
    int p1, p2, idx, jdx, kdx, ldx, index, position;
    Complex element;

    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_left_block = tensor_lattice->get_num_left_block();
    num_right_block = tensor_lattice->get_num_right_block();
    
    zero_flag = new bool** [left_bond];
    for(int o=0;o<left_bond;++o)
    {
        zero_flag[o] = new bool* [num_left_block];
        for(int j=0;j<num_left_block;++j)
        {
            zero_flag[o][j] = new bool[num_left_block];
            for(int i=0;i<num_left_block;++i)
                zero_flag[o][j][i] = false;
        }
    }

    // for each right index(block), find all possible left indiced(blocks)
    num_same_index = new int[num_right_block];
    position_same_index = new int* [num_right_block];
    for(int i=0;i<num_right_block;++i)
        tensor_lattice->ket_tensor_->FindMatrixBlock(num_same_index[i], position_same_index[i], 1, i);

    // if zero_flag[o1][j][i] is true, 
    // then right_contraction_tensor[o1]->find_matrix_block(j, i) is nonezero
    for(int o2=0;o2<right_bond;++o2)
    for(int i=0;i<right_contraction_tensor_[o2]->get_num_block();++i)
    {
        ldx = right_contraction_tensor_[o2]->left_index_[i];
        kdx = right_contraction_tensor_[o2]->right_index_[i];
        tmp_tensor[0] = right_contraction_tensor_[o2]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int l=0;l<num_same_index[ldx];++l)
        {
            jdx = tensor_lattice->ket_tensor_->left_index_[position_same_index[ldx][l]];
            p2 = tensor_lattice->physics_index_[jdx][ldx];
            if(tensor_operator->RightCheckZero(o2, p2) == false) continue; //zero
            for(int k=0;k<num_same_index[kdx];++k)
            {
                idx = tensor_lattice->ket_tensor_->left_index_[position_same_index[kdx][k]];
                p1 = tensor_lattice->physics_index_[idx][kdx];
                for(int o1=0;o1<left_bond;++o1)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0) zero_flag[o1][jdx][idx] = true;
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
            right_index = nullptr;
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
        delete tensor_contraction->right_contraction_tensor_[o1];
        tensor_contraction->right_contraction_tensor_[o1] = new ComplexMatrixBlock(num_block, left_index, right_index);

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
        ldx = right_contraction_tensor_[o2]->left_index_[i];
        kdx = right_contraction_tensor_[o2]->right_index_[i];
        tmp_tensor[0] = right_contraction_tensor_[o2]->get_matrix_block(i);
        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int l=0;l<num_same_index[ldx];++l)
        {
            jdx = tensor_lattice->ket_tensor_->left_index_[position_same_index[ldx][l]];
            p2 = tensor_lattice->physics_index_[jdx][ldx];
            if(tensor_operator->RightCheckZero(o2, p2) == false) continue; //zero
            // K*R (a_1, b_1)*(b_1, b_0) = (a_1, b_0)
            tmp_tensor[1] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[ldx][l]];
            first_tensor = tmp_tensor[1]->MultiplyToMatrix(tmp_tensor[0]);
            
            for(int k=0;k<num_same_index[kdx];++k)
            {
                idx = tensor_lattice->ket_tensor_->left_index_[position_same_index[kdx][k]];
                p1 = tensor_lattice->physics_index_[idx][kdx];
                // K*R*B (a_1, b_0)*(a_0, b_0)' = (a_1, a_0)
                tmp_tensor[2] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[kdx][k]];
                bra_tensor = tmp_tensor[2]->HermitianConjugateMatrix();
                second_tensor = first_tensor->MultiplyToMatrix(bra_tensor);
                delete bra_tensor;
                
                for(int o1=0;o1<left_bond;++o1)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0) 
                    {
                        // K*R*B*W
                        position = tensor_contraction->right_contraction_tensor_[o1]->FindMatrixBlock(jdx, idx);
                        
                        // if matrix is nullptr, replace it
                        tensor_contraction->right_contraction_tensor_[o1]->
                                            AddToMatrixBlock(position, element, second_tensor);
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
void ComplexTensorContraction::ComputeEffectHamilton(ComplexTensorLattice* tensor_lattice, 
        ComplexTensorOperator* tensor_operator, Complex* hamilton)
{
    ComplexMatrix *inv_tensor, *result_tensor, *tmp_tensor[2];
    Complex element, *matrix_element;
    int left_bond, right_bond;
    int vector_dim, position[2], result_dim[2], part_dim[2];
    int p1, p2, idx, jdx, kdx, ldx;

    vector_dim = tensor_lattice->ComputeKetTensorDim();

    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
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
            
            // if exist (p1,i,j) and (p2,k,l)
            position[0] = tensor_lattice->ket_tensor_->FindMatrixBlock(idx, jdx);
            position[1] = tensor_lattice->ket_tensor_->FindMatrixBlock(kdx, ldx);

            if(position[0]!=-1 && position[1]!=-1)
            {
                p1 = tensor_lattice->physics_index_[idx][jdx];
                p2 = tensor_lattice->physics_index_[kdx][ldx];
                element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                if(element.real_!=0.0 || element.imag_!=0.0)
                {
                    inv_tensor = tmp_tensor[1]->TransposeMatrix();
                    result_tensor = tmp_tensor[0]->MatrixKronProduct(inv_tensor);
                    delete inv_tensor;

                    if(result_tensor != nullptr)
                    {
                        result_tensor->MultiplyToScalar(element);

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
void ComplexTensorContraction::MultiplyEffectHamilton(ComplexTensorLattice* tensor_lattice, 
        ComplexTensorOperator* tensor_operator, Complex* state)
{
    ComplexMatrix *first_tensor, *second_tensor, *tmp_tensor[3];
    Complex *matrix_element, element;
    int left_bond, right_bond, num_right_block;
    int *num_same_index, **position_same_index, part_dim, position, total_element_num, p1, p2, idx, jdx, ldx, kdx;


    left_bond = tensor_operator->get_left_bond();
    right_bond = tensor_operator->get_right_bond();
    num_right_block = tensor_lattice->get_num_right_block();

    // for each right index(block), find all possible left indiced(blocks)
    num_same_index = new int[num_right_block];
    position_same_index = new int* [num_right_block];
    for(int i=0;i<num_right_block;++i)
        tensor_lattice->ket_tensor_->FindMatrixBlock(num_same_index[i], position_same_index[i], 1, i);

    // compute right_contraction_tensor
    for(int o2=0;o2<right_bond;++o2)
    for(int i=0;i<right_contraction_tensor_[o2]->get_num_block();++i)
    {
        ldx = right_contraction_tensor_[o2]->left_index_[i];
        jdx = right_contraction_tensor_[o2]->right_index_[i];
        tmp_tensor[0] = right_contraction_tensor_[o2]->get_matrix_block(i);

        if(tmp_tensor[0]->get_row() == 0) continue;
        for(int l=0;l<num_same_index[ldx];++l)
        {
            kdx = tensor_lattice->ket_tensor_->left_index_[position_same_index[ldx][l]];
            p2 = tensor_lattice->physics_index_[kdx][ldx];
            if(tensor_operator->RightCheckZero(o2, p2) == false) continue; //zero
            // K*R (a_1, b_1)*(b_1, b_0) = (a_1, b_0)
            tmp_tensor[1] = tensor_lattice->ket_tensor_->matrix_block_[position_same_index[ldx][l]];
            first_tensor = tmp_tensor[1]->MultiplyToMatrix(tmp_tensor[0]);
            for(int j=0;j<num_same_index[jdx];++j)
            {
                idx = tensor_lattice->ket_tensor_->left_index_[position_same_index[jdx][j]];
                p1 = tensor_lattice->physics_index_[idx][jdx];
                part_dim = tensor_lattice->ComputePartKetTensorDim(position_same_index[jdx][j]);
                for(int o1=0;o1<left_bond;++o1)
                {
                    element = tensor_operator->tensor_operator_[o1][o2]->get_matrix_element(p1, p2);
                    if(element.real_!=0.0 || element.imag_!=0.0) 
                    {
                        position = left_contraction_tensor_[o1]->FindMatrixBlock(idx, kdx);
                        if(position != -1)
                        {
                            tmp_tensor[2] = left_contraction_tensor_[o1]->matrix_block_[position];
                            second_tensor = tmp_tensor[2]->MultiplyToMatrix(first_tensor);
                            second_tensor->MultiplyToScalar(element);
                            matrix_element = second_tensor->get_matrix_element();
                            total_element_num = second_tensor->get_total_element_num();
                            for(int a=0;a<total_element_num;++a) state[part_dim+a] += matrix_element[a];
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
