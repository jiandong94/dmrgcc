#include "dmrg/realdmrg/real_tensor_lanczos.h"

RealTensorLanczos::RealTensorLanczos(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian, 
        RealTensorNetwork* network)
{
    space_ = space;
    hamiltonian_ = hamiltonian;
    network_ = network;
}

RealTensorLanczos::~RealTensorLanczos()
{}

void RealTensorLanczos::LanczosMethod(int num_iter, int site, double& result_value)
{
    RealTensorLattice* tensor_lattice;
    RealTensorOperator* tensor_operator;
    RealTensorContraction* tensor_contraction;
    double *eigenvector, *eigenvalue, *diagonal_element, *subdiagonal_element;
    double **psi, *h_psi, lanczos_precision, zero_tolerance;
    double *result_vector, old_energy, new_energy;
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
        eigenvector = new double[psi_dim*psi_dim];
        eigenvalue = new double[psi_dim];
        for(int i=0;i<psi_dim;++i) for(int j=0;j<psi_dim;++j)
        {
            eigenvector[i*psi_dim+j] = 0.0;
        }
        
        tensor_contraction->ComputeEffectHamilton(tensor_lattice, tensor_operator, eigenvector);
        
        RealSymMatrixDiag(eigenvector, eigenvalue, psi_dim);
        
        result_value = eigenvalue[0];
        result_vector = new double[psi_dim];
        for(int i=0;i<psi_dim;++i)
            result_vector[i] = eigenvector[i*psi_dim];
        tensor_lattice->VectorizeTensorLattice(false, result_vector);
        
        delete[] eigenvector;
        delete[] eigenvalue;
        delete[] result_vector;
    }
    else  // lanczos
    {
        eigenvector = new double[lanczos_iter*lanczos_iter];
        eigenvalue = new double[lanczos_iter];

        diagonal_element = new double[lanczos_iter];
        subdiagonal_element = new double[lanczos_iter];
        
        psi = new double* [lanczos_iter];
        for(int i=0;i<lanczos_iter;++i)
            psi[i] = nullptr;
        psi[0] = new double[psi_dim];
        h_psi = new double[psi_dim];
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
            subdiagonal_element[i] = sqrt(subdiagonal_element[i]);

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
            RealSymMatrixDiag(eigenvector, eigenvalue, lanczos_dim);
            new_energy = eigenvalue[0];
            
            if(subdiagonal_element[i]<zero_tolerance || 
               fabs(new_energy-old_energy)<lanczos_precision ||
               lanczos_dim == lanczos_iter)
            {
                result_value = eigenvalue[0];
                result_vector = new double[psi_dim];
                for(int j=0;j<psi_dim;++j)
                {
                    result_vector[j] = 0.0;
                    for(int k=0;k<lanczos_dim;++k)
                    {
                        result_vector[j] += psi[k][j]*eigenvector[k*lanczos_dim];
                    }
                }
                tensor_lattice->VectorizeTensorLattice(false, result_vector);
                double norm = 0;
                for(int i=0;i<psi_dim;++i)
                {
                    norm += result_vector[i]*result_vector[i];
                }
                if(fabs(norm-1.0) > 1E-6)
                {
                    cout << "false" << endl;
                }
                
                delete[] result_vector;
                break;
            }
            else
            {
                old_energy = new_energy;
            }

            if(lanczos_dim < lanczos_iter)
            {
                psi[i+1] = new double[psi_dim];
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

void RealTensorLanczos::VectorMultiply(double* vector1, double* vector2, int vector_dim, double &result)
{
    result = 0.0;
    for(int i=0;i<vector_dim;++i)
    {
        result += vector1[i]*vector2[i];
    }
}

void RealTensorLanczos::VectorSubtraction(double* vector1, double* vector2, int vector_dim, double factor2)
{
    for(int i=0;i<vector_dim;++i)
        vector1[i] -= factor2*vector2[i];
}

// \beta_n = \alpha_n - <\alpha_n, \beta_1>\beta_1 - 
//                ... - <\alpha_n, \beta_{n-1}>\beta_{n-1}
void RealTensorLanczos::GramSchmidtMethod(double** vector, int vector_dim, int num_vector)
{
    double result;
    for(int i=0;i<num_vector;++i)
    {
        VectorMultiply(vector[i], vector[num_vector], vector_dim, result);
        VectorSubtraction(vector[num_vector], vector[i], vector_dim, result);
    }

    VectorMultiply(vector[num_vector], vector[num_vector], vector_dim, result);
    result = sqrt(result);
    for(int i=0;i<vector_dim;++i)
        vector[num_vector][i] /= result;
}
