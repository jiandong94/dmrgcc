#ifndef DMRGCC_DMRG_COMPLEX_TENSOR_OPERATOR_H_
#define DMRGCC_DMRG_COMPLEX_TENSOR_OPERATOR_H_

#include "tensor/complex_matrix_block.h"

class ComplexTensorOperator
{
    friend class ComplexTensorContraction;
    protected:

    int physics_dim_;

    int left_bond_;
    int right_bond_;

    ComplexMatrix*** tensor_operator_;

    public:

    //
    //
    ComplexTensorOperator();

    //
    //
    ComplexTensorOperator(int physics_dim_);

    //
    //
    ComplexTensorOperator(ComplexMatrix* tmp_tensor, Complex coefficient);

    //
    //
    ~ComplexTensorOperator();

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
    void ExpanTensorOperator(ComplexMatrix** basic_operator, int leigh, int expan_operator_index, 
            Complex expan_coefficient);
    
    // deparallelisation algorithm
    //
    void LeftParallelTensorOperator(int &num_unparallel, int* &position_unparallel, 
            ComplexMatrix* &transfer_tensor);

    // deparallelisation algorithm
    //
    void RightParallelTensorOperator(int &num_unparallel, int* &position_unparallel, 
            ComplexMatrix* &transfer_tensor);

    // deparallelisation algorithm
    //
    void LeftMergeTensorOperator(int num_unparallel, ComplexMatrix* transfer_tensor);
    
    // deparallelisation algorithm
    //
    void RightMergeTensorOperator(int num_unparallel, ComplexMatrix* transfer_tensor);
    
    //
    //
    bool LeftCheckZero(int left_bond_index, int physics_index);
    
    //
    //
    bool RightCheckZero(int right_bond_index, int physics_index);

    //
    //
    bool MixedCheckZero(int left_bond_index, int physics_index);
};

inline void ComplexTensorOperator::DefineTensorOperator()
{
    tensor_operator_ = new ComplexMatrix** [left_bond_];
    for(int l=0;l<left_bond_;++l)
    {
        tensor_operator_[l] = new ComplexMatrix* [right_bond_];
        for(int r=0;r<right_bond_;++r)
            tensor_operator_[l][r] = new ComplexMatrix(physics_dim_, physics_dim_);
    }
}

inline void ComplexTensorOperator::ResetTensorOperator()
{
    for(int l=0;l<left_bond_;++l)
    {
        for(int r=0;r<right_bond_;++r) delete tensor_operator_[l][r];
        delete[] tensor_operator_[l];
    }
    delete[] tensor_operator_;
}
#endif // DMRGCC_DMRG_COMPLEX_TENSOR_OPERATOR_H_
