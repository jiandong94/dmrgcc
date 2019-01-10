#ifndef DMRGCC_DMRG_REAL_TENSOR_LANCZOS_H_
#define DMRGCC_DMRG_REAL_TENSOR_LANCZOS_H_

#include "dmrg/real_tensor_network.h"

class RealTensorLanczos
{
    protected:

    RealTensorSpace* space_;

    RealTensorHamiltonian* hamiltonian_;

    RealTensorNetwork* network_;

    bool disk_cache_;
    char cache_name_[512];

    int num_site_;
    int num_site_pp_;
    int num_site_mm_;

    int num_sweep_;
    int* num_iter_;
    int* max_block_;
    int* max_dim_;
    double* canonical_precision;
    double* noise_factor;

    public:

    RealTensorLanczos(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian, 
            int num_site, bool disk_cache);

};


#endif // DMRGCC_DMRG_REAL_TENSOR_LANCZOS_H_
