#include "tensor/real_matrix_block.h"


RealMatrixBlock::RealMatrixBlock()
{
    num_block_ = 0;
    left_index_ = nullptr;
    right_index_ = nullptr;

    matrix_block_ = nullptr;
}


RealMatrixBlock::RealMatrixBlock(int num_block)
{
   num_block_ = num_block;

   if(num_block_ > 0)
   {
       left_index_ = new int [num_block_];
       right_index_ = new int [num_block_];

       matrix_block_ = new RealMatrix* [num_block_];

       for(int i=0;i<num_block_;++i)
       {
           left_index_[i] = i;
           right_index_[i] = i;
           matrix_block_[i] = new RealMatrix();
       }
   }
   else
   {
       left_index_ = nullptr;
       right_index_ = nullptr;
       matrix_block_ = nullptr;
   }

}


RealMatrixBlock::RealMatrixBlock(int num_block, int* left_index, int* right_index)
{
   num_block_ = num_block;

   if(num_block_ > 0)
   {
       left_index_ = new int [num_block_];
       right_index_ = new int [num_block_];

       matrix_block_ = new RealMatrix* [num_block_];

       for(int i=0;i<num_block_;++i)
       {
           left_index_[i] = left_index[i];
           right_index_[i] = right_index[i];
           matrix_block_[i] = new RealMatrix();
       }
   }
   else
   {
       left_index_ = nullptr;
       right_index_ = nullptr;
       matrix_block_ = nullptr;
   }

}


RealMatrixBlock::RealMatrixBlock(RealMatrixBlock* tmp_matrix_block)
{
   num_block_ = tmp_matrix_block->num_block_;

   if(num_block_ > 0)
   {
       left_index_ = new int [num_block_];
       right_index_ = new int [num_block_];

       matrix_block_ = new RealMatrix* [num_block_];

       for(int i=0;i<num_block_;++i)
       {
           left_index_[i] = tmp_matrix_block->left_index_[i];
           right_index_[i] = tmp_matrix_block->right_index_[i];
           matrix_block_[i] = tmp_matrix_block->matrix_block_[i];
       }
   }
   else
   {
       left_index_ = nullptr;
       right_index_ = nullptr;
       matrix_block_ = nullptr;
   }

}

RealMatrixBlock::~RealMatrixBlock()
{
    if(matrix_block_ != nullptr)
    {
        delete[] left_index_;
        delete[] right_index_;

        for(int i=0;i<num_block_;++i) delete matrix_block_[i];

        delete[] matrix_block_;
    }
}

int RealMatrixBlock::get_num_block()
{
    return num_block_;   
}

int* RealMatrixBlock::get_left_index()
{
    return left_index_;
}

int* RealMatrixBlock::get_right_index()
{
    return right_index_;
}

RealMatrix* RealMatrixBlock::get_matrix_block(int position)
{
    return matrix_block_[position];
}

void RealMatrixBlock::set_matrix_block(int position, RealMatrix* tmp_matrix)
{
    if(position < 0 || position > num_block_)
    {
        cout << "index is out of range" << endl;
        exit(-1);
    }
    else
    {
        matrix_block_[position] = new RealMatrix(tmp_matrix);
    }
}

int RealMatrixBlock::ComputeMatrixBlockDim()
{
    int matrix_block_dim = 0;
    if(matrix_block_ != nullptr)
    {
        for(int i=0;i<num_block_;++i)
            matrix_block_dim += matrix_block_[i]->get_total_element_num();
    }
    return matrix_block_dim;
}

int RealMatrixBlock::ComputePartMatrixBlockDim(int position)
{
    int part_matrix_block_dim = 0;
    for(int i=0;i<position;++i)
    {
        part_matrix_block_dim += matrix_block_[i]->get_total_element_num();
    }
    return part_matrix_block_dim;
}

void RealMatrixBlock::NormalizeMatrixBlock()
{
    double result = 0.0;
    for(int i=0;i<num_block_;++i) result += matrix_block_[i]->SumSquareMatrix();

    result = sqrt(1.0/result);
    for(int i=0;i<num_block_;++i) matrix_block_[i]->MultiplyToScalar(result);
}

void RealMatrixBlock::VectorizeMatrixBlock(bool direction, double* state)
{
    double* matrix_element;
    int total_element_num, part_matrix_block_dim;
    for(int i=0;i<num_block_;++i)
    {
        matrix_element = matrix_block_[i]->get_matrix_element();
        total_element_num = matrix_block_[i]->get_total_element_num();
        part_matrix_block_dim = ComputePartMatrixBlockDim(i);
        if(direction == true) //to state
        {    
            for(int j=0;j<total_element_num;++j)
                state[part_matrix_block_dim+j] = matrix_element[j];
        }
        else if(direction == false) //to block
        {
            for(int j=0;j<total_element_num;++j)
                matrix_element[j] = state[part_matrix_block_dim+j];
        }
    }
    
}

void RealMatrixBlock::PrintMatrixBlock()
{
    if(matrix_block_ == nullptr)
    {
        cout << "matrix_block_ is nullptr!" << endl;
    }
    else
    {
        cout << "==============================" << endl;
        cout << "Print MatrixBlock: " << endl;
        cout << "Number of blocks: " << num_block_ << endl;
        for(int i=0;i<num_block_;++i)
        {
            cout << "Left index= " << left_index_[i] << ", " <<
                "Right index= " << right_index_[i] << endl;
            matrix_block_[i]->PrintMatrix();
        }
    }
}

void RealMatrixBlock::WriteMatrixBlock(const char* matrix_block_name)
{
    ofstream matrix_block_file;

    matrix_block_file.open(matrix_block_name, ios::binary|ios::out);

    WriteMatrixBlock(matrix_block_file);

    matrix_block_file.close();
}

void RealMatrixBlock::WriteMatrixBlock(ofstream &matrix_block_file)
{
    matrix_block_file.write((char*) &num_block_, sizeof(int));
    if(matrix_block_ != nullptr)
    {
       for(int i=0;i<num_block_;++i)
       {
            matrix_block_file.write((char*) &left_index_[i], sizeof(int));
            matrix_block_file.write((char*) &right_index_[i], sizeof(int));
            matrix_block_[i]->WriteMatrix(matrix_block_file);
       }
    }
}


void RealMatrixBlock::ReadMatrixBlock(const char* matrix_block_name)
{
    ifstream matrix_block_file;

    matrix_block_file.open(matrix_block_name, ios::binary|ios::out);
    
    ReadMatrixBlock(matrix_block_file);

    matrix_block_file.close();
}


void RealMatrixBlock::ReadMatrixBlock(ifstream &matrix_block_file)
{
    if(matrix_block_ != nullptr)
    {
        delete[] left_index_;
        delete[] right_index_;
        
        for(int i=0;i<num_block_;++i) delete matrix_block_[i];

        delete[] matrix_block_;
    }

    matrix_block_file.read((char*) &num_block_, sizeof(int));
    if(num_block_ > 0)
    {
        left_index_ = new int[num_block_];
        right_index_ = new int[num_block_];
        matrix_block_ = new RealMatrix*[num_block_];

        for(int i=0;i<num_block_;++i)
        {
            matrix_block_file.read((char*) &left_index_[i], sizeof(int));
            matrix_block_file.read((char*) &right_index_[i], sizeof(int));

            matrix_block_[i] = new RealMatrix();
            matrix_block_[i]->ReadMatrix(matrix_block_file);
        }
    }

}



void RealMatrixBlock::ResetMatrixBlock()
{
    for(int i=0;i<num_block_;++i) matrix_block_[i]->ResetMatrix();

}

void RealMatrixBlock::AddToMatrixBlock(int position, RealMatrix* tmp_matrix)
{
    if(matrix_block_ != nullptr)
    {
        matrix_block_[position]->AddToMatrix(tmp_matrix);
    }
}

void RealMatrixBlock::MultiplyToScalar(double scalar)
{
    for(int i=0;i<num_block_;++i)
    {
        if(matrix_block_[i] == nullptr)
        {
            cout << "nullptr*scalar in RealMatrixBlock::MultiplyScalar" << endl;
            exit(-1);
        }
        matrix_block_[i]->MultiplyToScalar(scalar);
    }
}

int RealMatrixBlock::FindMatrixBlock(int left, int right)
{
    int position;

    position = -1;
    for(int i=0;i<num_block_;++i)
    {
        if(left_index_[i] == left && right_index_[i] == right)
        {
            position = i;
            break;
        }

    }
    return position;
}


void RealMatrixBlock::FindMatrixBlock(int &num_same_index, int* &position_same_index,
                                      int leigh, int target_index)
{
    num_same_index = 0;
    if(leigh == 0)
    {
        for(int i=0;i<num_block_;++i)
        {
            if(left_index_[i]==target_index) num_same_index++;
        }
    }
    else if(leigh == 1)
    {
        for(int i=0;i<num_block_;++i)
        {
            if(right_index_[i]==target_index) num_same_index++;
        }
    }

    int p = 0;
    position_same_index = nullptr;
    if(num_same_index > 0)
    {
        position_same_index = new int[num_same_index];
        if(leigh == 0)
        {
            for(int i=0;i<num_block_;++i)
            {
                if(left_index_[i]==target_index)
                {
                    position_same_index[p] = i;
                    p++;
                }
            }
                
        }
        else if(leigh == 1)
        {
            for(int i=0;i<num_block_;++i)
            {
                if(right_index_[i]==target_index)
                {
                    position_same_index[p] = i;
                    p++;
                }
            }
        }
    }
}

































