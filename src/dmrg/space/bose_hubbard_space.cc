#include "dmrg/space/bose_hubbard_space.h"

BoseHubbardSpace::BoseHubbardSpace(int num_site_x, int num_site_y, int num_boson, int physics_dim)
{
    // derive class member
    num_site_x_ = num_site_x;
    num_site_y_ = num_site_y;
    num_boson_ = num_boson;
    physics_dim_ = physics_dim;
    
    // base class member
    disk_cache_ = false;
    num_site_ = num_site_x_*num_site_y_;
    num_site_pp_ = num_site_x_*num_site_y_+1;
    num_site_mm_ = num_site_x_*num_site_y_-1;
    tensor_lattice_ = new RealTensorLattice* [num_site_];
    for(int i=0;i<num_site_;++i)
        tensor_lattice_[i] = new RealTensorLattice(physics_dim_);
    DefineQuantumTable();
}

BoseHubbardSpace::~BoseHubbardSpace()
{

}

void BoseHubbardSpace::DefineQuantumTable()
{
    double mean_quantum;
    int min_from_left, max_from_left, min_from_right, max_from_right, min_final, max_final;
    num_quantum_ = 1;
    num_table_ = new int[num_site_pp_];
    quantum_table_ = new int** [num_site_pp_];

    for(int i=0;i<num_site_pp_;++i)
    {
        if(i == 0)
        {
            min_final = 0;
            max_final = 0;
            mean_quantum = 0;
        }
        else
        {
            min_from_left = 0;
            max_from_left = (physics_dim_-1)*i;
            min_from_right = num_boson_-(num_site_-i)*(physics_dim_-1);
            max_from_right = num_boson_;
            min_final = max(min_from_left, min_from_right);
            max_final = min(max_from_left, max_from_right);

            mean_quantum = num_boson_/num_site_*i;
        }
        num_table_[i] = max_final - min_final+1;
        quantum_table_[i] = new int* [num_table_[i]];
        for(int j=0;j<num_table_[i];++j)
            quantum_table_[i][j] = new int[num_quantum_];
        for(int j=min_final;j<max_final+1;++j)
            quantum_table_[i][j-min_final][0] = j;
        ReorderQuantumTable(num_quantum_, num_table_[i], quantum_table_[i], &mean_quantum);
    }

}

void BoseHubbardSpace::MergeQuantumTable(int* merge_quantum_table, int* operator_quantum_table, 
        int* space_quantum_table)
{
    merge_quantum_table[0] = operator_quantum_table[0] + space_quantum_table[0];
}

int BoseHubbardSpace::CheckQuantumTable(int site, int* left_table, int* right_table)
{
    int info = -1;

    for(int i=0;i<physics_dim_;++i)
    {
        if(left_table[0]+i == right_table[0])
        {
            info = i;
            break;
        }
    }
    return info;
}
