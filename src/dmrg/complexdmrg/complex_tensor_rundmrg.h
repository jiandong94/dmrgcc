#ifndef DMRGCC_DMRG_COMPLEX_TENSOR_RUNDMRG_H_
#define DMRGCC_DMRG_COMPLEX_TENSOR_RUNDMRG_H_

#include "dmrg/complexdmrg/complex_tensor_lanczos.h"
#include "util/input.h"
class ComplexTensorRundmrg
{
    protected:

    ComplexTensorSpace* space_;

    ComplexTensorHamiltonian* hamiltonian_;

    ComplexTensorNetwork* network_;

    bool disk_cache_;
    char cache_name_[512];
    bool cache_record_;
    char record_name_[512];
    bool cache_resume_;
    char resume_name_[512];
    bool record_process_;
    char process_name_[512];


    int num_site_;
    int num_site_pp_;
    int num_site_mm_;

    int num_sweep_;
    int* max_dim_;
    int* max_block_;
    int* num_iter_;
    double* canonical_precision_;
    double* noise_factor_;

    public:

    //
    //
    ComplexTensorRundmrg(ComplexTensorSpace* space, ComplexTensorHamiltonian* hamiltonian,
        InputGroup& input);

    //
    //
    ~ComplexTensorRundmrg();

    //
    //
    void Sweep(InputGroup &table);

    //
    //
    void Run();

    // right canonical tensor space and compute
    // all right tensor
    void Initialize();
};

#endif // DMRGCC_DMRG_COMPLEX_TENSOR_RUNDMRG_H_
