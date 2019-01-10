#include "dmrg/real_tensor_network.h"

RealTensorNetwork::RealTensorNetwork(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian)
{
    int left_bond, right_bond;

    space_ = space;
    hamiltonian_ = hamiltonian;

    disk_cache_ = space->get_disk_cache();
    strcpy(cache_name_, space->get_cache_name());

    num_site_ = hamiltonian->num_site_;
    num_site_pp_ = hamiltonian->num_site_pp_;
    num_site_mm_ = hamiltonian->num_site_mm_;

    tensor_contraction_ = new RealTensorContraction* [num_site_pp_];
    for(int i=0;i<num_site_pp_;++i)
    {
        if(i<num_site_)
        {
            left_bond = hamiltonian->tensor_hamiltonian_[i]->get_left_bond();
            right_bond = hamiltonian->tensor_hamiltonian_[i]->get_right_bond();
        }
        else
        {
            left_bond = hamiltonian->tensor_hamiltonian_[num_site_mm_]->get_right_bond();
            right_bond = hamiltonian->tensor_hamiltonian_[0]->get_left_bond();
        }

        tensor_contraction_[i] = new RealTensorContraction(left_bond, right_bond);
    }
    
}

RealTensorNetwork::~RealTensorNetwork()
{
    for(int i=0;i<num_site_pp_;++i)
        delete tensor_contraction_[i];
    delete[] tensor_contraction_;
}

RealTensorContraction* RealTensorNetwork::get_tensor_contraction(int site)
{
    return tensor_contraction_[site];
}

void RealTensorNetwork::DefineTensorNetwork()
{
    tensor_contraction_[0]->DefineTensorContraction(0);
    tensor_contraction_[num_site_mm_]->DefineTensorContraction(1);

    if(disk_cache_ == true)
    {
        RecordTensorNetwork(0, 0);
        ResetTensorNetwork(0, 0);

        RecordTensorNetwork(1, num_site_mm_);
        ResetTensorNetwork(1, num_site_mm_);
    }
}

void RealTensorNetwork::ResetTensorNetwork(int leigh, int site)
{
    tensor_contraction_[site]->ResetTensorContraction(leigh);
}

void RealTensorNetwork::RecordTensorNetwork(int leigh, int site)
{
    char tensor_name[512];

    sprintf(tensor_name, "%s/NetworkDiskCacheFile/TensorNetwork_leigh_%d_site_%d.dat", 
                          cache_name_, leigh, site);

    tensor_contraction_[site]->WriteTensorContraction(leigh, tensor_name);
}

void RealTensorNetwork::ResumeTensorNetwork(int leigh, int site)
{
    char tensor_name[512];

    sprintf(tensor_name, "%s/NetworkDiskCacheFile/TensorNetwork_leigh_%d_site_%d.dat", 
                          cache_name_, leigh, site);

    tensor_contraction_[site]->ReadTensorContraction(leigh, tensor_name);
}

void RealTensorNetwork::RemoveTensorNetwork(int leigh, int site)
{
    char command[512];
    int info;

    if(disk_cache_ == true)
    {
        sprintf(command, " rm -rf %s/NetworkDiskCacheFile/TensorNetwork_leigh_%d_site_%d.dat", 
                          cache_name_, leigh, site);
        info = system(command);
    }
}

void RealTensorNetwork::ExpanTensorNetwork(int leigh, int site)
{
    int operator_num_table, **operator_quantum_table, **mapping_table;
    
    operator_num_table = hamiltonian_->num_table_[site+(leigh+1)%2];
    operator_quantum_table = hamiltonian_->quantum_table_[site+(leigh+1)%2];

    space_->ComputeExpanTensorLattice((leigh+1)%2, site, operator_num_table, operator_quantum_table, 
                                mapping_table);
    if(leigh==0 && site<num_site_mm_)
    {
        tensor_contraction_[site]->LeftExpanTensorContraction(space_->tensor_lattice_[site], 
               hamiltonian_->tensor_hamiltonian_[site], space_->expan_tensor_lattice_, 
               mapping_table, space_->noise_factor_);
    }
    else if(leigh==1 && site>0)
    {
        tensor_contraction_[site]->RightExpanTensorContraction(space_->tensor_lattice_[site], 
               hamiltonian_->tensor_hamiltonian_[site], space_->expan_tensor_lattice_, 
               mapping_table, space_->noise_factor_);
    }

    for(int i=0;i<operator_num_table;++i)
        delete[] mapping_table[i];
    delete mapping_table;
}

void RealTensorNetwork::ComputeTensorNetwork(int leigh, int site)
{
    if(leigh==0) 
    {
        tensor_contraction_[site]->LeftComputeTensorContraction(space_->tensor_lattice_[site], 
                hamiltonian_->tensor_hamiltonian_[site], tensor_contraction_[site+1]);
    }
    else if(leigh==1) 
    {
        tensor_contraction_[site]->RightComputeTensorContraction(space_->tensor_lattice_[site], 
                hamiltonian_->tensor_hamiltonian_[site], tensor_contraction_[(site-1+num_site_pp_)%num_site_pp_]);
    }
}

void RealTensorNetwork::ComputeTensorNetwork()
{
    DefineTensorNetwork();
    
    for(int i=0;i<num_site_;++i)
    {
        ComputeTensorNetwork(0, i);
        ResetTensorNetwork(0, i);
    }

}
