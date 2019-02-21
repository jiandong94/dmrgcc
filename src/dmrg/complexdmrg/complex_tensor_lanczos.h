#ifndef DMRGCC_DMRG_COMPLEX_TENSOR_LANCZOS_H_
#define DMRGCC_DMRG_COMPLEX_TENSOR_LANCZOS_H_

#include "dmrg/complexdmrg/complex_tensor_network.h"

class ComplexTensorLanczos
{
    protected:

    ComplexTensorSpace* space_;

    ComplexTensorHamiltonian* hamiltonian_;

    ComplexTensorNetwork* network_;

    public:

    ComplexTensorLanczos(ComplexTensorSpace* space, ComplexTensorHamiltonian* hamiltonian, 
            ComplexTensorNetwork* network);

    ~ComplexTensorLanczos();

    void LanczosMethod(int num_times_sweep, int site, double& result);

    void VectorMultiply(Complex* vector1, Complex* vector2, int vector_dim, Complex &result);

    void VectorSubtraction(Complex* vector1, Complex* vector2, int vector_dim, Complex factor2);

    void GramSchmidtMethod(Complex **vector, int vector_dim, int num_vector);

};


#endif // DMRGCC_DMRG_COMPLEX_TENSOR_LANCZOS_H_
