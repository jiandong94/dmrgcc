#include "dmrg/complexdmrg/complex_tensor_lanczos.h"

ComplexTensorLanczos::ComplexTensorLanczos(ComplexTensorSpace* space, ComplexTensorHamiltonian* hamiltonian, 
        ComplexTensorNetwork* network)
{
    space_ = space;
    hamiltonian_ = hamiltonian;
    network_ = network;
}

ComplexTensorLanczos::~ComplexTensorLanczos()
{}

void ComplexTensorLanczos::LanczosMethod(int num_iter, int site, double& result_value)
{
    ComplexTensorLattice* tensor_lattice;
    ComplexTensorOperator* tensor_operator;
    ComplexTensorContraction* tensor_contraction;
    Complex *eigenvector, *result_vector, *diagonal_element, *subdiagonal_element;
    Complex **psi, *h_psi;
    double lanczos_precision, zero_tolerance, *eigenvalue, old_energy, new_energy;
    int psi_dim, lanczos_iter, lanczos_dim;

    tensor_lattice = space_->get_tensor_lattice(site);
    tensor_operator = hamiltonian_->get_tensor_hamiltonian(site);
    tensor_contraction = network_->get_tensor_contraction(site);
    
    lanczos_iter = num_iter;
    lanczos_precision = 1E-15;
    zero_tolerance = 1E-15;

    psi_dim = tensor_lattice->ComputeKetTensorDim();
    // exact diagonalize
    if(psi_dim < 100)
    {
        eigenvector = new Complex[psi_dim*psi_dim];
        eigenvalue = new double[psi_dim];
        for(int i=0;i<psi_dim;++i) for(int j=0;j<psi_dim;++j)
        {
            eigenvector[i*psi_dim+j] = 0.0;
        }
        
        tensor_contraction->ComputeEffectHamilton(tensor_lattice, tensor_operator, eigenvector);
        ComplexSymMatrixDiag(eigenvector, eigenvalue, psi_dim);
        
        result_value = eigenvalue[0];
        result_vector = new Complex[psi_dim];
        for(int i=0;i<psi_dim;++i)
            result_vector[i] = eigenvector[i*psi_dim];
        tensor_lattice->VectorizeTensorLattice(false, result_vector);
        
        delete[] eigenvector;
        delete[] eigenvalue;
        delete[] result_vector;
    }
    else  // lanczos
    {
        eigenvector = new Complex[lanczos_iter*lanczos_iter];
        eigenvalue = new double[lanczos_iter];

        diagonal_element = new Complex[lanczos_iter];
        subdiagonal_element = new Complex[lanczos_iter];
        
        psi = new Complex* [lanczos_iter];
        for(int i=0;i<lanczos_iter;++i)
            psi[i] = nullptr;
        psi[0] = new Complex[psi_dim];
        h_psi = new Complex[psi_dim];
        for(int i=0;i<psi_dim;++i)
        {
            psi[0][i] = 0.0;
            h_psi[i] = 0.0;
        }
        tensor_lattice->VectorizeTensorLattice(true, psi[0]);

        old_energy = 10000.0;

        for(int i=0;i<lanczos_iter;++i)
        {
            //tensor_lattice->PrintTensorLattice();
            tensor_contraction->MultiplyEffectHamilton(tensor_lattice, tensor_operator, h_psi);
            // d = <psi|H|psi>
            VectorMultiply(psi[i], h_psi, psi_dim, diagonal_element[i]);
            // h_psi = H|psi> - d|psi>
            VectorSubtraction(h_psi, psi[i], psi_dim, diagonal_element[i]);
            // f = <h_psi|h_psi>^{-1}
            VectorMultiply(h_psi, h_psi, psi_dim, subdiagonal_element[i]);
            subdiagonal_element[i].real_ = sqrt(subdiagonal_element[i].real_);

            // tridiagonal matrix
            lanczos_dim = i+1;
            for(int j=0;j<lanczos_dim;++j)
            {
                for(int k=0;k<lanczos_dim;++k)
                {
                    eigenvector[k+j*lanczos_dim] = 0.0;
                }
                eigenvector[j+j*lanczos_dim] = diagonal_element[j];
                if(j > 0)
                {
                    eigenvector[j+(j-1)*lanczos_dim] = subdiagonal_element[j-1];
                    eigenvector[(j-1)+j*lanczos_dim] = subdiagonal_element[j-1];
                }
            }
            // diagonalize the tridiagonal matrix
            ComplexSymMatrixDiag(eigenvector, eigenvalue, lanczos_dim);
            new_energy = eigenvalue[0];
            
            if(subdiagonal_element[i]<zero_tolerance || 
               fabs(new_energy-old_energy)<lanczos_precision ||
               lanczos_dim == lanczos_iter)
            {
                result_value = eigenvalue[0];
                result_vector = new Complex[psi_dim];
                for(int j=0;j<psi_dim;++j)
                {
                    result_vector[j] = 0.0;
                    for(int k=0;k<lanczos_dim;++k)
                    {
                        result_vector[j] += psi[k][j]*eigenvector[k*lanczos_dim];
                    }
                }
                tensor_lattice->VectorizeTensorLattice(false, result_vector);
                
                delete[] result_vector;
                break;
            }
            else
            {
                old_energy = new_energy;
            }

            if(lanczos_dim < lanczos_iter)
            {
                psi[i+1] = new Complex[psi_dim];
                for(int j=0;j<psi_dim;++j)
                {
                    psi[i+1][j] = h_psi[j]/subdiagonal_element[i];
                    // ?h_psi[j] = diagonal_element[i]*psi[i][j];
                    h_psi[j] = -subdiagonal_element[i]*psi[i][j];
                }

                GramSchmidtMethod(psi, psi_dim, lanczos_dim);

                tensor_lattice->VectorizeTensorLattice(false, psi[i+1]);
            }
        }

        delete[] eigenvector;
        delete[] eigenvalue;
        delete[] diagonal_element;
        delete[] subdiagonal_element;

        for(int i=0;i<lanczos_iter;++i)
            delete[] psi[i];
        delete[] psi;
        delete[] h_psi;
    }
}

void ComplexTensorLanczos::VectorMultiply(Complex* vector1, Complex* vector2, int vector_dim, Complex &result)
{
    result = 0.0;
    for(int i=0;i<vector_dim;++i)
        result += vector1[i]*vector2[i];
}

void ComplexTensorLanczos::VectorSubtraction(Complex* vector1, Complex* vector2, int vector_dim, Complex factor2)
{
    for(int i=0;i<vector_dim;++i)
        vector1[i] -= factor2*vector2[i];
}

// \beta_n = \alpha_n - <\alpha_n, \beta_1>\beta_1 - 
//                ... - <\alpha_n, \beta_{n-1}>\beta_{n-1}
void ComplexTensorLanczos::GramSchmidtMethod(Complex** vector, int vector_dim, int num_vector)
{
    Complex result;
    for(int i=0;i<num_vector;++i)
    {
        VectorMultiply(vector[i], vector[num_vector], vector_dim, result);
        VectorSubtraction(vector[num_vector], vector[i], vector_dim, result);
    }

    VectorMultiply(vector[num_vector], vector[num_vector], vector_dim, result);
    result.real_ = sqrt(result.real_);
    for(int i=0;i<vector_dim;++i)
        vector[num_vector][i] /= result.real_;
}
