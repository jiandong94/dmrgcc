#include "dmrg/real_tensor_operator.h"

RealTensorOperator::RealTensorOperator()
{
    physics_dim_ = 0;
    
    left_bond_ = 0;
    right_bond_ = 0;

    tensor_operator_ = nullptr;
}

RealTensorOperator::RealTensorOperator(int physics_dim)
{
    physics_dim_ = physics_dim;

    left_bond_ = 0;
    right_bond_ = 0;

    tensor_operator_ = nullptr;
}

RealTensorOperator::RealTensorOperator(RealMatrix* tmp_tensor, double coefficient)
{
    double element;
    physics_dim_ = tmp_tensor->get_row();

    left_bond_ = 1;
    right_bond_ = 1;

    DefineTensorOperator();
    for(int i=0;i<physics_dim_;++i) for(int j=0;j<physics_dim_;++j)
    {
        element = tmp_tensor->get_matrix_element(i, j);
        tensor_operator_[0][0]->set_matrix_element(i, j, coefficient*element);
    }
}

RealTensorOperator::~RealTensorOperator()
{
    ResetTensorOperator();
}

int RealTensorOperator::get_physics_dim()
{
    return physics_dim_;
}

int RealTensorOperator::get_left_bond()
{
    return left_bond_;
}

int RealTensorOperator::get_right_bond()
{
    return right_bond_;
}

void RealTensorOperator::PrintTensorOperator()
{
    cout << "==============================" << endl;
    cout << "RealTensorOperator: " << endl;
    cout << "Physics dimension: " << physics_dim_ << endl;
    cout << "Left bond dimension: " << left_bond_ <<
           " Right bond dimension: " << right_bond_ << endl;
    cout << endl;
    cout << "Operator: " << endl;
    if(tensor_operator_ == nullptr)
    {
        cout << "operator_tensor is nullptr!" << endl;
    }
    else
    {
        for(int i=0;i<left_bond_;++i) 
        {
            for(int j=0;j<right_bond_;++j)
            {
                cout << "left bond = " << i << ", right bond = " << j << endl;
                tensor_operator_[i][j]->PrintMatrix();
            }
        }
    }
}

void RealTensorOperator::WriteTensorOperator(const char* tensor_operator_name)
{
    ofstream tensor_operator_file;

    tensor_operator_file.open(tensor_operator_name, ios::binary|ios::out);  

    WriteTensorOperator(tensor_operator_file);

    tensor_operator_file.close();
}

void RealTensorOperator::WriteTensorOperator(ofstream &tensor_operator_file)
{
    tensor_operator_file.write((char*) &physics_dim_, sizeof(int));

    tensor_operator_file.write((char*) &left_bond_, sizeof(int));
    tensor_operator_file.write((char*) &right_bond_, sizeof(int));

    for(int i=0;i<left_bond_;++i) for(int j=0;j<right_bond_;++j)
        tensor_operator_[i][j]->WriteMatrix(tensor_operator_file);
}


void RealTensorOperator::ReadTensorOperator(const char* tensor_operator_name)
{
    ifstream tensor_operator_file;

    tensor_operator_file.open(tensor_operator_name, ios::binary|ios::out);  

    ReadTensorOperator(tensor_operator_file);

    tensor_operator_file.close();
}


void RealTensorOperator::ReadTensorOperator(ifstream &tensor_operator_file)
{
    ResetTensorOperator();
    
    tensor_operator_file.read((char*) &physics_dim_, sizeof(int));

    tensor_operator_file.read((char*) &left_bond_, sizeof(int));
    tensor_operator_file.read((char*) &right_bond_, sizeof(int));

    tensor_operator_ = new RealMatrix** [left_bond_];
    for(int i=0;i<left_bond_;++i) 
    {
        tensor_operator_[i] = new RealMatrix* [right_bond_];
        for(int j=0;j<right_bond_;++j)
        {
            tensor_operator_[i][j] = new RealMatrix(physics_dim_, physics_dim_);
            tensor_operator_[i][j]->ReadMatrix(tensor_operator_file);
        }
    }
}

void RealTensorOperator::ExpanTensorOperator(RealMatrix** basic_operator, int leigh, 
        int expan_operator_index, double expan_coefficient)
{
    RealMatrix*** expan_tensor_operator;
    int expan_left_bond, expan_right_bond, edge[2], expan_index[2];

    edge[0] = 1;
    edge[1] = 1;

    if(leigh == 0) edge[0] = 0;
    if(leigh == 1) edge[1] = 0;

    if(left_bond_==0 && right_bond_==0)
    {
        expan_left_bond = 1;
        expan_right_bond = 1;
    }
    else
    {
        expan_left_bond = left_bond_ + edge[0];
        expan_right_bond = right_bond_ + edge[1];
    }
    
    // expan tensor operator
    // first site: right_bond_++
    // last site: left_bond_++
    // middle sites: left_bond_++ right_bond_++
    expan_tensor_operator = new RealMatrix** [expan_left_bond];
    for(int i=0;i<expan_left_bond;++i)
    {
        expan_tensor_operator[i] = new RealMatrix* [expan_right_bond];
        for(int j=0;j<expan_right_bond;++j)
        {
            expan_tensor_operator[i][j] = new RealMatrix(physics_dim_, physics_dim_);
        }
    }

    // first  site: tensor_operator[0][right_bond_] = new operator
    // last   site: tensor_operator[left_bond_][0] = new operator
    // middle site: tensor_operator[left_bond_][right_bond_] = new operator
    for(int i=0;i<physics_dim_;++i) for(int j=0;j<physics_dim_;++j)
    {
        for(int l=0;l<left_bond_;++l) for(int r=0;r<right_bond_;++r)
        {
            expan_tensor_operator[l][r]->set_matrix_element(i, j, 
                    tensor_operator_[l][r]->get_matrix_element(i, j));
        }
        expan_index[0] = edge[0]*left_bond_;
        expan_index[1] = edge[1]*right_bond_;
        expan_tensor_operator[expan_index[0]][expan_index[1]]->set_matrix_element(i, j, 
        expan_coefficient*basic_operator[expan_operator_index]->get_matrix_element(i, j));
    }
    ResetTensorOperator();
    left_bond_ = expan_left_bond;
    right_bond_ = expan_right_bond;
    tensor_operator_ = expan_tensor_operator;
}

void RealTensorOperator::LeftParallelTensorOperator(int &num_unparallel, int* &position_unparallel, 
        RealMatrix* &transfer_tensor)
{
    RealMatrix* reshape_tensor;
    double element;
    int reshape_row, reshape_column, reshape_position;

    reshape_row = left_bond_ * physics_dim_ * physics_dim_;
    reshape_column = right_bond_;
    reshape_tensor = new RealMatrix(reshape_row, reshape_column);

    // reshape_row (p2,p1,l)
    for(int l=0;l<left_bond_;++l) for(int r=0;r<right_bond_;++r)
    for(int p1=0;p1<physics_dim_;++p1) for(int p2=0;p2<physics_dim_;++p2)
    {   
        // row first
        reshape_position = p2+(p1+l*physics_dim_)*physics_dim_;
        element = tensor_operator_[l][r]->get_matrix_element(p1, p2);
        reshape_tensor->set_matrix_element(reshape_position, r, element);
    }

    reshape_tensor->ParallelMatrix(0, num_unparallel, position_unparallel, transfer_tensor);

    ResetTensorOperator();
    right_bond_ = num_unparallel;
    DefineTensorOperator();

    for(int l=0;l<left_bond_;++l) for(int r=0;r<right_bond_;++r)
    for(int p1=0;p1<physics_dim_;++p1) for(int p2=0;p2<physics_dim_;++p2)
    {   
        // row first
        reshape_position = p2+(p1+l*physics_dim_)*physics_dim_;
        element = reshape_tensor->get_matrix_element(reshape_position, position_unparallel[r]);
        tensor_operator_[l][r]->set_matrix_element(p1, p2, element);
    }

    delete reshape_tensor;
}

void RealTensorOperator::RightParallelTensorOperator(int &num_unparallel, int* &position_unparallel, 
        RealMatrix* &transfer_tensor)
{
    RealMatrix* reshape_tensor;
    double element;
    int reshape_row, reshape_column, reshape_position;

    reshape_row = left_bond_;
    reshape_column = right_bond_ * physics_dim_ * physics_dim_;
    reshape_tensor = new RealMatrix(reshape_row, reshape_column);

    // reshape_row (p2,p1,r)
    for(int l=0;l<left_bond_;++l) for(int r=0;r<right_bond_;++r)
    for(int p1=0;p1<physics_dim_;++p1) for(int p2=0;p2<physics_dim_;++p2)
    {   
        // row first
        reshape_position = p2+(p1+l*physics_dim_)*physics_dim_;
        element = tensor_operator_[l][r]->get_matrix_element(p1, p2);
        reshape_tensor->set_matrix_element(l, reshape_position, element);
    }

    reshape_tensor->ParallelMatrix(1, num_unparallel, position_unparallel, transfer_tensor);

    ResetTensorOperator();
    left_bond_ = num_unparallel;
    DefineTensorOperator();

    for(int l=0;l<left_bond_;++l) for(int r=0;r<right_bond_;++r)
    for(int p1=0;p1<physics_dim_;++p1) for(int p2=0;p2<physics_dim_;++p2)
    {   
        // row first
        reshape_position = p2+(p1+l*physics_dim_)*physics_dim_;
        element = reshape_tensor->get_matrix_element(position_unparallel[l], reshape_position);
        tensor_operator_[l][r]->set_matrix_element(p1, p2, element);
    }

    delete reshape_tensor;
}

void RealTensorOperator::LeftMergeTensorOperator(int num_unparallel, RealMatrix* transfer_tensor)
{
    RealMatrix ***result_tensor_operator;
    double element, factor[2];
    int merge_dim;

    merge_dim = transfer_tensor->get_column();
    if(merge_dim != left_bond_)
    {
        error("merge dimension is not match in LeftMergeTensorOperator");
    }

    result_tensor_operator = new RealMatrix** [num_unparallel];
    for(int l=0;l<num_unparallel;++l)
    {
        result_tensor_operator[l] = new RealMatrix*[right_bond_];
        for(int r=0;r<right_bond_;++r)
            result_tensor_operator[l][r] = new RealMatrix(physics_dim_, physics_dim_);
    }

    //      |
    // --M--O--
    //      |
    for(int l=0;l<num_unparallel;++l) for(int r=0;r<right_bond_;++r)
    for(int p1=0;p1<physics_dim_;++p1) for(int p2=0;p2<physics_dim_;++p2)
    {
        element = 0.0;
        
        for(int m=0;m<merge_dim;++m)
        {
            factor[0] = transfer_tensor->get_matrix_element(l, m);
            factor[1] = tensor_operator_[m][r]->get_matrix_element(p1, p2);
            element += factor[0]*factor[1];
        }

        result_tensor_operator[l][r]->set_matrix_element(p1, p2, element);
    }
    
    ResetTensorOperator();
    left_bond_ = num_unparallel;
    tensor_operator_ = result_tensor_operator;
}

void RealTensorOperator::RightMergeTensorOperator(int num_unparallel, RealMatrix* transfer_tensor)
{
    RealMatrix ***result_tensor_operator;
    double element, factor[2];
    int merge_dim;

    merge_dim = transfer_tensor->get_row();
    if(merge_dim != right_bond_)
    {
        error("merge dimension is not match in RightMergeTensorOperator");
    }

    result_tensor_operator = new RealMatrix** [left_bond_];
    for(int l=0;l<left_bond_;++l)
    {
        result_tensor_operator[l] = new RealMatrix*[num_unparallel];
        for(int r=0;r<num_unparallel;++r)
            result_tensor_operator[l][r] = new RealMatrix(physics_dim_, physics_dim_);
    }

    //   |
    // --O--M--
    //   |
    for(int l=0;l<left_bond_;++l) for(int r=0;r<num_unparallel;++r)
    for(int p1=0;p1<physics_dim_;++p1) for(int p2=0;p2<physics_dim_;++p2)
    {
        element = 0.0;
        
        for(int m=0;m<merge_dim;++m)
        {
            factor[0] = transfer_tensor->get_matrix_element(m, r);
            factor[1] = tensor_operator_[l][m]->get_matrix_element(p1, p2);
            element += factor[0]*factor[1];
        }

        result_tensor_operator[l][r]->set_matrix_element(p1, p2, element);
    }
    
    ResetTensorOperator();
    right_bond_ = num_unparallel;
    tensor_operator_ = result_tensor_operator;
}

bool RealTensorOperator::LeftCheckZero(int l, int p1)
{
    bool zero_flag = false;
    for(int r=0;r<right_bond_;++r)
        for(int p2=0;p2<physics_dim_;++p2)
        {
            if(fabs(tensor_operator_[l][r]->get_matrix_element(p1, p2)) > 1E-15)
            {
                zero_flag = true;
                break;
            }
        }
    
    return zero_flag;
}

bool RealTensorOperator::RightCheckZero(int r, int p2)
{
    bool zero_flag = false;
    for(int l=0;l<left_bond_;++l)
        for(int p1=0;p1<physics_dim_;++p1)
        {
            if(fabs(tensor_operator_[l][r]->get_matrix_element(p1, p2)) > 1E-15)
            {
                zero_flag = true;
                break;
            }
        }
    
    return zero_flag;
}

bool RealTensorOperator::MixedCheckZero(int l, int p2)
{
    bool zero_flag = false;
    for(int r=0;r<right_bond_;++r)
        for(int p1=0;p1<physics_dim_;++p1)
        {
            if(fabs(tensor_operator_[l][r]->get_matrix_element(p1, p2)) > 1E-15)
            {
                zero_flag = true;
                break;
            }
        }
    
    return zero_flag;
}
