#ifndef DMRGCC_DMRG_COMPLEX_TENSOR_NETWORK_H_
#define DMRGCC_DMRG_COMPLEX_TENSOR_NETWORK_H_

#include "dmrg/complexdmrg/complex_tensor_contraction.h"
#include "dmrg/complexdmrg/complex_tensor_space.h"
#include "dmrg/complexdmrg/complex_tensor_hamiltonian.h"

class ComplexTensorNetwork
{
    protected:

    ComplexTensorSpace* space_;

    ComplexTensorHamiltonian* hamiltonian_;

    bool disk_cache_;
    char cache_name_[1024];

    int num_site_;
    int num_site_pp_;
    int num_site_mm_;

    ComplexTensorContraction** tensor_contraction_;

    public:

    //
    //
    ComplexTensorNetwork(ComplexTensorSpace* space, ComplexTensorHamiltonian* hamiltonian);

    //
    //
    ~ComplexTensorNetwork();

    //
    //
    ComplexTensorContraction* get_tensor_contraction(int site);

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

#endif //#ifndef DMRGCC_DMRG_COMPLEX_TENSOR_NETWORK_H_
