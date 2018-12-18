#ifndef DMRGCC_DMRG_REAL_TENSOR_LATTICE_H_
#define DMRGCC_DMRG_REAL_TENSOR_LATTICE_H_

#include "tensor/real_matrix_block.h"

class RealTensorLattice
{
    protected:
    // physics dimension
    int physics_dim_;

    // left side
    int num_left_block_;
    int* left_block_;
    int* left_dim_;

    // right side
    int num_right_block_;
    int* right_block_;
    int* right_dim_;

    // physics index
    int** physics_index_;

    // braket tensor
    RealMatrixBlock* ket_tensor_;

    // match dimension
    int* match_dim_;

    // canonical tensor
    RealMatrixBlock* canonical_tensor_;

    public:
    // constructor
    //
    RealTensorLattice();

    // constructor
    //
    RealTensorLattice(int physics_dim);

    // constructor
    //
    RealTensorLattice(RealTensorLattice* tmp_tensor_lattice);

    // destructor
    //
    ~RealTensorLattice();

    //
    //
    int get_physics_dim();

    //
    //
    int get_num_left_block();

    //
    //
    int* get_left_block();

    //
    //
    int* get_left_dim();
    
    //
    //
    int get_num_right_block();

    //
    //
    int* get_right_block();

    //
    //
    int* get_right_dim();

    // 
    //
    int** get_physics_dim();

    //
    //
    int ComputeLatticeDim(int leigh);

    //
    //
    int ComputeKetTensorDim();

    //
    //
    void PrintTensorLattice();

    //
    //
    void WriteTensorLattice(char* tensor_lattice_name);

    //
    //
    void WriteTensorLattice(ofstream &tensor_lattice_file);

    //
    //
    void ReadTensorLattice(char* tensor_lattice_name);

    //
    //
    void ReadTensorLattice(ifstream &tensor_lattice_file);


















};





#endif // DMRGCC_DMRG_REAL_TENSOR_LATTICE_H_
