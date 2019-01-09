#ifndef DMRGCC_SPACE_BOSE_HUBBARD_SPACE_H_
#define DMRGCC_SPACE_BOSE_HUBBARD_SPACE_H_

#include "dmrg/real_tensor_space.h"

class BoseHubbardSpace : public RealTensorSpace
{
    protected:

    int num_site_x_;
    int num_site_y_;

    int num_boson_;

    int num_level_;

    public:

    // constructor
    //
    BoseHubbardSpace(int num_site_x, int num_site_y, int num_boson, int num_level);

    //
    //
    ~BoseHubbardSpace();

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



#endif // DMRGCC_SPACE_BOSE_HUBBARD_SPACE_H_
