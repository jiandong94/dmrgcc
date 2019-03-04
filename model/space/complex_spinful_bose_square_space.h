#ifndef DMRGCC_SPACE_COMPLEX_SPINFUL_BOSE_SQUARE_SPACE_H_
#define DMRGCC_SPACE_COMPLEX_SPINFUL_BOSE_SQUARE_SPACE_H_

#include "dmrg/complexdmrg/complex_tensor_space.h"

class ComplexSpinfulBoseSquareSpace : public ComplexTensorSpace
{
    protected:

    int num_site_x_;
    int num_site_y_;

    int num_boson_[2];

    int physics_dim_;

    public:

    // constructor
    //
    ComplexSpinfulBoseSquareSpace(int num_site_x, int num_site_y, int* num_boson, int physics_dim);

    //
    //
    ~ComplexSpinfulBoseSquareSpace();

    protected:

    //
    //
    void DefineQuantumTable();

    //
    //
    void MergeQuantumTable(int* merge_quantum_table, int* operator_quantum_table, 
            int* space_quantum_table);

    //
    //
    int CheckQuantumTable(int site, int* left_table, int* right_table);

};



#endif // DMRGCC_SPACE_COMPLEX_SPINFUL_BOSE_SQUARE_SPACE_H_
