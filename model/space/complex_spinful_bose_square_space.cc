#include "model/space/complex_spinful_bose_square_space.h"

ComplexSpinfulBoseSquareSpace::ComplexSpinfulBoseSquareSpace(int num_site_x, int num_site_y, int* num_boson, int physics_dim)
{
    // derive class member
    num_site_x_ = num_site_x;
    num_site_y_ = num_site_y;
    num_boson_[0] = num_boson[0];
    num_boson_[1] = num_boson[1];
    physics_dim_ = physics_dim;
    
    // base class member
    num_site_ = 2*num_site_x_*num_site_y_;
    num_site_pp_ = 2*num_site_x_*num_site_y_+1;
    num_site_mm_ = 2*num_site_x_*num_site_y_-1;
    
    tensor_lattice_ = new ComplexTensorLattice* [num_site_];
    for(int i=0;i<num_site_;++i)
        tensor_lattice_[i] = new ComplexTensorLattice(physics_dim_);
    DefineQuantumTable();
}

ComplexSpinfulBoseSquareSpace::~ComplexSpinfulBoseSquareSpace()
{

}

void ComplexSpinfulBoseSquareSpace::DefineQuantumTable()
{
    double mean_quantum[2];
    int min_from_left[2], max_from_left[2], min_from_right[2], max_from_right[2],
        min_final[2], max_final[2], actual_position[2];
    int actual_site, tmp_table_index;
    num_quantum_ = 2;
    num_table_ = new int[num_site_pp_];
    quantum_table_ = new int** [num_site_pp_];

    for(int i=0;i<num_site_pp_;++i)
    {
        if(i == 0)
        {
            min_final[0] = 0;
            min_final[1] = 0;
            max_final[0] = 0;
            max_final[1] = 0;
            mean_quantum[0] = 0;
            mean_quantum[1] = 0;
        }
        else
        {
            // actual position for up: (i+1)/2
            // actual position for down: i/2
            actual_position[0] = (i+1)/2;
            actual_position[1] = i/2;
            actual_site = num_site_/2;

            min_from_left[0] = 0;
            max_from_left[0] = (physics_dim_-1)*actual_position[0];
            min_from_left[1] = 0;
            max_from_left[1] = (physics_dim_-1)*actual_position[1];
            
            min_from_right[0] = num_boson_[0]-(actual_site-actual_position[0])*(physics_dim_-1);
            max_from_right[0] = num_boson_[0];
            min_from_right[1] = num_boson_[1]-(actual_site-actual_position[1])*(physics_dim_-1);
            max_from_right[1] = num_boson_[1];
            
            for(int t=0;t<2;++t)
            {
                min_final[t] = max(min_from_left[t], min_from_right[t]);
                max_final[t] = min(max_from_left[t], max_from_right[t]);
            }
            
            mean_quantum[0] = (double)num_boson_[0]/actual_site*actual_position[0];
            mean_quantum[1] = (double)num_boson_[1]/actual_site*actual_position[1];
        }

        num_table_[i] = 1;
        for(int j=0;j<2;++j)
           num_table_[i] *= max_final[j] - min_final[j]+1;

        quantum_table_[i] = new int* [num_table_[i]];
        for(int j=0;j<num_table_[i];++j)
            quantum_table_[i][j] = new int[num_quantum_];
        
        for(int t0=min_final[0];t0<max_final[0]+1;++t0)
        for(int t1=min_final[1];t1<max_final[1]+1;++t1)
        {
            tmp_table_index = (t0-min_final[0]) + (t1-min_final[1])*(max_final[0]-min_final[0]+1);
            quantum_table_[i][tmp_table_index][0] = t0;
            quantum_table_[i][tmp_table_index][1] = t1;
        }
        ReorderQuantumTable(num_quantum_, num_table_[i], quantum_table_[i], mean_quantum);
    }

}

void ComplexSpinfulBoseSquareSpace::MergeQuantumTable(int* merge_quantum_table, int* operator_quantum_table, 
        int* space_quantum_table)
{
    merge_quantum_table[0] = operator_quantum_table[0] + space_quantum_table[0];
    merge_quantum_table[1] = operator_quantum_table[1] + space_quantum_table[1];
}

int ComplexSpinfulBoseSquareSpace::CheckQuantumTable(int site, int* left_table, int* right_table)
{
    int info = -1;
    int shift, physics_dim[2];
    
    shift = site%2;

    for(int i=0;i<physics_dim_;++i)
    {
        physics_dim[0] = 0;
        physics_dim[1] = 0;
        physics_dim[shift] = i;

        if(left_table[0]+physics_dim[0] == right_table[0])
        if(left_table[1]+physics_dim[1] == right_table[1])
        {
            info = i;
            break;
        }
    }
    return info;
}
