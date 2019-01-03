#ifndef DMRGCC_DMRG_REAL_TENSOR_OPERATOR_H_
#define DMRGCC_DMRG_REAL_TENSOR_OPERATOR_H_

#include "itensor/matrix_block.h"

class RealTensorOperator
{
    protected:

    int physics_dim_;

    int left_bond_;
    int right_bond;

    RealTensorMatrix*** tensor_operator_;

    public:

    //
    //
    RealTensorOperator();

    //
    //
    RealTensorOperator(int physics_dim_);

    //
    //
    RealTensorOperator(RealMatrix* tmp_tensor, double coefficient);

    //
    //
    ~RealTensorOperator();

    //
    //
    int get_physics_dim();

    //
    //
    int get_left_bond();

    //
    //
    int get_right_bond();

    //
    //
    void PrintTensorOperator();

    //
    //
    void WriteTensorOperator(const char* tensor_operator_name);

    //
    //
    void WriteTensorOperator(ofstream &tensor_operator_file);

    //
    //
    void ReadTensorOperator(const char* tensor_operator_name);

    //
    //
    void ReadTensorOperator(ifstream &tensor_operator_file);

    //
    //
    void DefineTensorOperator();

    //
    //
    void ResetTensorOperator();

    //
    //
    void ExpanTensorOperator(RealMatrix** basic_operator, int leigh, int expan_operator_index, 
            double expan_coefficient);

    //
    //
    bool LeftCheckZero(int left_bond_index, int physics_index);
    
    //
    //
    bool RightCheckZero(int right_bond_index, int physics_index);
};

inline void RealTensorOperator::DefineTensorOperator()
{
    tensor_operator_ = new RealMatrix** [left_bond_];
    for(int l=0;l<left_bond_;++l)
    {
        tensor_operator_[l] = new RealMatrix* [right_bond_];
        for(int r=0;r<right_bond_;++r)
            tensor_operator_[l][r] = RealMatrix(physics_dim_);

    }
}

inline void RealTensorOperator::ResetTensorOperator()
{
    for(int l=0;l<left_bond_;++l)
    {
        for(int r=0;r<right_bond_;++r) delete tensor_operator_[l][r];
        delete[] tensor_operator[l];
    }
    delete[] tensor_operator;
}
#endif // DMRGCC_DMRG_REAL_TENSOR_OPERATOR_H_
