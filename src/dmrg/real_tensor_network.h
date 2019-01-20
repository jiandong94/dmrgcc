#ifndef DMRGCC_DMRG_REAL_TENSOR_NETWORK_H_
#define DMRGCC_DMRG_REAL_TENSOR_NETWORK_H_

#include "dmrg/real_tensor_contraction.h"
#include "dmrg/real_tensor_space.h"
#include "dmrg/real_tensor_hamiltonian.h"

class RealTensorNetwork
{
    protected:

    RealTensorSpace* space_;

    RealTensorHamiltonian* hamiltonian_;

    bool disk_cache_;
    char cache_name_[512];

    int num_site_;
    int num_site_pp_;
    int num_site_mm_;

    RealTensorContraction** tensor_contraction_;

    public:

    //
    //
    RealTensorNetwork(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian);

    //
    //
    ~RealTensorNetwork();

    //
    //
    RealTensorContraction* get_tensor_contraction(int site);

    //
    //
    void DefineTensorNetwork();

    //
    //
    void ResetTensorNetwork(int leigh, int site);

    //
    //
    void PrintTensorNetwork();

    //
    //
    void RecordTensorNetwork(int leigh, int site);

    //
    //
    void ResumeTensorNetwork(int leigh, int site);

    //
    //
    void RemoveTensorNetwork(int leigh, int site);

    //
    //
    void ExpanTensorNetwork(int leigh, int site);

    //
    //
    void ComputeTensorNetwork(int leigh, int site);

    //
    //
    void ComputeTensorNetwork();
};

#endif //#ifndef DMRGCC_DMRG_REAL_TENSOR_NETWORK_H_
