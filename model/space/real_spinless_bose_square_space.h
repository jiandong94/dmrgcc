#ifndef DMRGCC_SPACE_REAL_SPINLESS_BOSE_SQUARE_SPACE_H_
#define DMRGCC_SPACE_REAL_SPINLESS_BOSE_SQUARE_SPACE_H_

#include "dmrg/realdmrg/real_tensor_space.h"

class RealSpinlessBoseSquareSpace : public RealTensorSpace
{
    protected:

    int num_site_x_;
    int num_site_y_;

    int num_boson_;

    int physics_dim_;

    public:

    // constructor
    //
    RealSpinlessBoseSquareSpace(int num_site_x, int num_site_y, int num_boson, int physics_dim);

    //
    //
    ~RealSpinlessBoseSquareSpace();

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



#endif // DMRGCC_SPACE_REAL_SPINLESS_BOSE_SQUARE_SPACE_H_
