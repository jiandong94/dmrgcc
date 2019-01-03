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
    if(tensor_tensor_ == nullptr)
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

    tensor_operator_file.open(operator_name, ios::binary|ios::out);  

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
        tensor_operator_[l] = new RealMatrix* [right_bond_];
        for(int j=0;j<right_bond_;++j)
        {
            tensor_operator_[i][j] = new RealMatrix(physics_dim_, physics_dim_);
            tensor_operator_[i][j]->ReadMatrix(tensor_operator_file);
        }
    }
}

void RealTensorOperator::ExpanTensorOperator(RealMatrix** basic_operator, int leigh, 
        int expan_operator_index, double expan_coefficient);
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
        expan_right_bond = right_bond_ = edge[1];
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
                    tensor_operator[l][r]->get_matrix_element(i, j));
        }
        expan_index[0] = edge[0]*left_bond_;
        expan_index[1] = edge[1]*right_bond_;
        expan_tensor_operator[expan_index[0]][expan_index[1]]->set_matrix_element(i, j, 
        expan_coefficient*basic_tensor[expan_operator_index]->get_matrix_element(i, j));
    }
    ResetTensorOperator();
    left_bond_ = expan_left_bond;
    right_bond_ = expan_right_bond;
    tensor_operator_ = expan_tensor_operator;
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
