#ifndef DMRGCC_DMRG_REAL_TENSOR_LANCZOS_H_
#define DMRGCC_DMRG_REAL_TENSOR_LANCZOS_H_

#include "dmrg/real_tensor_network.h"

class RealTensorLanczos
{
    protected:

    RealTensorSpace* space_;

    RealTensorHamiltonian* hamiltonian_;

    RealTensorNetwork* network_;

    public:

    RealTensorLanczos(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian, 
            RealTensorNetwork* network);

    ~RealTensorLanczos();

    void LanczosMethod(int num_times_sweep, int site, double& result);

    void VectorMultiply(double* vector1, double* vector2, int vector_dim, double &result);

    void VectorSubtraction(double* vector1, double* vector2, int vector_dim, double factor2);

    void GramSchmidtMethod(double **vector, int vector_dim, int num_vector);

};


#endif // DMRGCC_DMRG_REAL_TENSOR_LANCZOS_H_
