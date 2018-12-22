#ifndef DMRGCC_DMRG_REAL_TENSOR_LATTICE_H_
#define DMRGCC_DMRG_REAL_TENSOR_LATTICE_H_

#include "tensor/real_matrix_block.h"

class RealTensorLattice
{
    protected:
    // physics dimension
    int physics_dim_;

    // left side
    // num_left_block_ donates the number of the left
    // QuantumNumber, the block number in left_block_ is different
    // from each other and array left_dim_ stores dimensions for 
    // each block number.
    int num_left_block_;
    int* left_block_;
    int* left_dim_;

    // right side
    // similar to left side.
    int num_right_block_;
    int* right_block_;
    int* right_dim_;

    // physics index
    // for every left block number and right block number ,there
    // is a physics number stored in table physics_index_.
    int** physics_index_;

    // braket tensor
    // ket_tensor_ is a series of matrices with all possible left block 
    // number and right block number.
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
    int** get_physics_index();

    //
    //
    int ComputeLatticeDim(int leigh);

    //
    //
    int ComputeKetTensorDim();

    //
    //
    int ComputePartKetTensorDim(int position);

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

    //
    //
    void DefineTensorLattice(int num_left_block, int num_right_block, 
         int* left_block, int* right_block, int* left_dim, int* right_dim);

    //
    //
    void DefineTensorLattice(int num_block, int* left_index, int* right_index, 
         int* physics_index, int row, int column);


    // combine two tensor lattice
    // if leigh=0, two tensor lattice must have the same
    // num_right_block_ and right_block_, the new tensor
    // lattice contains the left_block_ combined by the two, 
    // if they have the same right block index, then add up
    // the left dimension.
    // eg: for leigh=0 
    // original tensor lattice: num_left_block_ = 3
    //                          left_block_ = [1,3,5]
    //                          left_dim_ = [2,3,4]
    // temp tensor lattice    : num_left_block_ = 4
    //                          left_block_ = [1,3,4,6]
    //                          left_dim_ = [1,3,4,6]
    // new tensor lattice     : num_left_block = 5
    //                          left_block_ = [1,3,4,5,6]
    //                          left_dim_ = [3,6,4,4,6]
    //void CombineTensorLattice(int leigh, int &num_leigh_block, int* &leigh_block, 
    //                          int* &leigh_dim, RealTensorLattice* tmp_tensor_lattice);
    //

    //
    //
    void ResetTensorLattice();

    //
    //
    void NormalizeTensorLattice();

    //
    //
    void VectorizeTensorLattice(bool direction, double* &state);

    //
    //
    void ComputeTruncateDim(int max_dim, double canonical_precision, int num_singular_block, 
         int* singular_dim, double** singular_value, int* truncate_dim);
    //
    //
    void LeftCanonicalTensorLattice(int max_dim, double canonical_precision);

    //
    //
    void LeftCanonicalTensorLattice(int* &singular_dim, int** &singular_value);

    //
    //
    void RightCanonicalTensorLattice(int max_dim, double canonical_precision);

    //
    //
    void RightCanonicalTensorLattice(int* &singular_dim, int** &singular_value);

    //
    //
    void MixCanonicalTensorLattice(int max_dim, double canonical_precision);

    //
    //
    void mixCanonicalTensorLattice(int* &singular_dim, int** &singular_value);












};





#endif // DMRGCC_DMRG_REAL_TENSOR_LATTICE_H_
