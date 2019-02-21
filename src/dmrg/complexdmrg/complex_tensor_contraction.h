#ifndef DMRGCC_DMRG_COMPLEX_TENSOR_CONTRACTION_H_
#define DMRGCC_DMRG_COMPLEX_TENSOR_CONTRACTION_H_

#include "tensor/complex_matrix_block.h"
#include "dmrg/complexdmrg/complex_tensor_lattice.h"
#include "dmrg/complexdmrg/complex_tensor_operator.h"

class ComplexTensorContraction
{
    protected:

    int left_bond_;
    int right_bond_;

    ComplexMatrixBlock** left_contraction_tensor_;
    ComplexMatrixBlock** right_contraction_tensor_;

    public:

    //
    //
    ComplexTensorContraction(int left_bond, int right_bond);

    //
    //
    ~ComplexTensorContraction();

    //
    //
    int get_left_bond();

    //
    //
    int get_right_bond();

    //
    //
    ComplexMatrixBlock** get_left_contraction_tensor();

    //
    //
    ComplexMatrixBlock** get_right_contraction_tensor();

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
    void LeftExpanTensorContraction(ComplexTensorLattice* tensor_lattice, ComplexTensorOperator* tensor_operator, 
            ComplexTensorLattice* expan_tensor_lattice, int** mapping_table, double noise_factor);

    //
    //
    void RightExpanTensorContraction(ComplexTensorLattice* tensor_lattice, ComplexTensorOperator* tensor_operator, 
            ComplexTensorLattice* expan_tensor_lattice, int** mapping_table, double noise_factor);

    //
    //
    void LeftComputeTensorContraction(ComplexTensorLattice* tensor_lattice, ComplexTensorOperator* tensor_operator, 
            ComplexTensorContraction* tensor_contraction);

    //
    //
    void RightComputeTensorContraction(ComplexTensorLattice* tensor_lattice, ComplexTensorOperator* tensor_operator, 
            ComplexTensorContraction* tensor_contraction);

    //
    //
    void ComputeEffectHamilton(ComplexTensorLattice* tensor_lattice, ComplexTensorOperator* tensor_operator, 
            Complex* hamilton);

    //
    //
    void MultiplyEffectHamilton(ComplexTensorLattice* tensor_lattice, ComplexTensorOperator* tensor_operator, 
            Complex* state);
};

#endif //DMRGCC_DMRG_COMPLEX_TENSOR_VACTION_H_
