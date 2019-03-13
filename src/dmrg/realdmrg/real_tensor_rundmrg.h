#ifndef DMRGCC_DMRG_REAL_TENSOR_RUNDMRG_H_
#define DMRGCC_DMRG_REAL_TENSOR_RUNDMRG_H_

#include "dmrg/realdmrg/real_tensor_lanczos.h"
#include "util/input.h"
class RealTensorRundmrg
{
    protected:

    RealTensorSpace* space_;

    RealTensorHamiltonian* hamiltonian_;

    RealTensorNetwork* network_;

    bool disk_cache_;
    char cache_name_[1024];
    bool cache_record_;
    char record_name_[1024];
    bool cache_resume_;
    char resume_name_[1024];
    bool record_process_;
    char process_name_[1024];


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
    RealTensorRundmrg(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian,
        InputGroup& input);

    //
    //
    ~RealTensorRundmrg();

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

#endif // DMRGCC_DMRG_REAL_TENSOR_RUNDMRG_H_
