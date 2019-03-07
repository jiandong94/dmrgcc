#ifndef DMRGCC_DMRG_COMPLEX_TENSOR_LATTICE_H_
#define DMRGCC_DMRG_COMPLEX_TENSOR_LATTICE_H_

#include "tensor/complex_matrix_block.h"

class ComplexTensorLattice
{
    friend class ComplexTensorContraction;
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
    ComplexMatrixBlock* ket_tensor_;

    // match dimension
    int* match_dim_;

    // canonical tensor
    ComplexMatrixBlock* canonical_tensor_;

    public:
    // constructor
    //
    ComplexTensorLattice();

    // constructor
    //
    ComplexTensorLattice(int physics_dim);

    // constructor
    //
    ComplexTensorLattice(ComplexTensorLattice* tmp_tensor_lattice);

    // destructor
    //
    ~ComplexTensorLattice();

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
    ComplexMatrixBlock* get_ket_tensor();

    //
    //
    ComplexMatrixBlock* get_canonical_tensor();

    //
    //
    int* get_match_dim();
    
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
    void WriteTensorLattice(const char* tensor_lattice_name);

    //
    //
    void WriteTensorLattice(ofstream &tensor_lattice_file);

    //
    //
    void ReadTensorLattice(const char* tensor_lattice_name);

    //
    //
    void ReadTensorLattice(ifstream &tensor_lattice_file);

    //
    //
    void DefineTensorLattice(const int num_left_block, const int num_right_block, 
         const int* left_block, const int* right_block, const int* left_dim, const int* right_dim);

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
    void CombineTensorLattice(int leigh, int &num_leigh_block, int* &leigh_block, 
                              int* &leigh_dim, ComplexTensorLattice* expan_tensor_lattice);
    

    //
    //
    void ResetTensorLattice();

    //
    //
    void NormalizeTensorLattice();

    //
    //
    void VectorizeTensorLattice(bool direction, Complex* &state);

    //
    //
    void ComputeTruncateDim(int max_dim, double canonical_precision, int num_singular_block, 
         int* singular_dim, double** singular_value, int* truncate_dim);
    //
    //
    void LeftCanonicalTensorLattice(int max_dim, double canonical_precision);

    //
    //
    void LeftCanonicalTensorLattice(int* &singular_dim, double** &singular_value);

    //
    //
    void RightCanonicalTensorLattice(int max_dim, double canonical_precision);


    void MixCanonicalTensorLattice(int* &singular_dim, int** &singular_value);

    //
    //
    void LeftMergeTensorLattice(ComplexTensorLattice* tmp_tensor_lattice);

    //
    //
    void RightMergeTensorLattice(ComplexTensorLattice* tmp_tensor_lattice);

    //
    //
    void LeftExpanTensorLattice(const ComplexTensorLattice* tmp_tensor_lattice, 
            const ComplexTensorLattice* expan_tensor_lattice);

    //
    //
    void RightExpanTensorLattice(ComplexTensorLattice* tmp_tensor_lattice, 
            ComplexTensorLattice* expan_tensor_lattice);
};

#endif // DMRGCC_DMRG_COMPLEX_TENSOR_LATTICE_H_
