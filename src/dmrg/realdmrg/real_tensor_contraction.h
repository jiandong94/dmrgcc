#ifndef DMRGCC_DMRG_REAL_TENSOR_CONTRACTION_H_
#define DMRGCC_DMRG_REAL_TENSOR_CONTRACTION_H_

#include "tensor/real_matrix_block.h"
#include "dmrg/realdmrg/real_tensor_lattice.h"
#include "dmrg/realdmrg/real_tensor_operator.h"

class RealTensorContraction
{
    protected:

    int left_bond_;
    int right_bond_;

    RealMatrixBlock** left_contraction_tensor_;
    RealMatrixBlock** right_contraction_tensor_;

    public:

    //
    //
    RealTensorContraction(int left_bond, int right_bond);

    //
    //
    ~RealTensorContraction();

    //
    //
    int get_left_bond();

    //
    //
    int get_right_bond();

    //
    //
    RealMatrixBlock** get_left_contraction_tensor();

    //
    //
    RealMatrixBlock** get_right_contraction_tensor();

    //
    //
    void PrintTensorContraction();

    //
    //
    void WriteTensorContraction(int leigh, char* tensor_contraction_name);

    //
    //
    void ReadTensorContraction(int leigh, char* tensor_contraction_name);

    //
    //
    void DefineTensorContraction(int leigh);

    //
    //
    void ResetTensorContraction(int leigh);

    //
    //
    void LeftExpanTensorContraction(RealTensorLattice* tensor_lattice, RealTensorOperator* tensor_operator, 
            RealTensorLattice* expan_tensor_lattice, int** mapping_table, double noise_factor);

    //
    //
    void RightExpanTensorContraction(RealTensorLattice* tensor_lattice, RealTensorOperator* tensor_operator, 
            RealTensorLattice* expan_tensor_lattice, int** mapping_table, double noise_factor);

    //
    //
    void LeftComputeTensorContraction(RealTensorLattice* tensor_lattice, RealTensorOperator* tensor_operator, 
            RealTensorContraction* tensor_contraction);

    //
    //
    void RightComputeTensorContraction(RealTensorLattice* tensor_lattice, RealTensorOperator* tensor_operator, 
            RealTensorContraction* tensor_contraction);

    //
    //
    void ComputeEffectHamilton(RealTensorLattice* tensor_lattice, RealTensorOperator* tensor_operator, 
            double* hamilton);

    //
    //
    void MultiplyEffectHamilton(RealTensorLattice* tensor_lattice, RealTensorOperator* tensor_operator, 
            double* state);
};

















#endif //DMRGCC_DMRG_REAL_TENSOR_VACTION_H_
