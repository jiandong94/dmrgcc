#ifndef DMRGCC_DMRG_REAL_TENSOR_RUNDMRG_H_
#define DMRGCC_DMRG_REAL_TENSOR_RUNDMRG_H_

#include "dmrg/real_tensor_network.h"
#include "util/input.h"
class RealTensorRundmrg
{
    protected:

    RealTensorSpace* space_;

    RealTensorHamiltonian* hamiltonian_;

    RealTensorNetwork* network_;

    bool disk_cache_;
    string cache_name_;

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
            bool disk_cache, string cache_name, int num_sweep, InputGroup table);

    //
    //
    ~RealTensorRundmrg();

    //
    //
    void Sweep(InputGroup &table);
};

#endif // DMRGCC_DMRG_REAL_TENSOR_RUNDMRG_H_
