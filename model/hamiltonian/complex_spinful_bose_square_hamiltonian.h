#ifndef DMRGCC_DMRG_COMPLEX_SPINFUL_BOSE_SQUARE_HAMILTONIAN_H_
#define DMRGCC_DMRG_COMPLEX_SPINFUL_BOSE_SQUARE_HAMILTONIAN_H_

#include "dmrg/complexdmrg/complex_tensor_hamiltonian.h"
class ComplexSpinfulBoseSquareHamiltonian : public ComplexTensorHamiltonian
{
    protected:
    
    int num_site_x_;
    int num_site_y_;
    int physics_dim_;

    double flux_[2];
    double hop_x_[2];
    double hop_y_[2];
    double inter_value_[3];
    double chemical_value_[2];
    
    bool period_x_;
    bool period_y_;

    public:

    //
    //
    ComplexSpinfulBoseSquareHamiltonian(int num_site_x, int num_site_y, int physics_dim, double* flux, 
            double* hop_x, double* hop_y, double* inter_value, double* chemical_value, 
            bool period_x=false, bool period_y=false);

    //
    //
    ~ComplexSpinfulBoseSquareHamiltonian();

    //
    //
    void DefineBasicTensor();

    //
    //
    void DefineHamiltonian1(int p1, int p2, Complex coefficient);

    //
    //
    void DefineHamiltonian2(int p1, int p2, Complex coefficient);
    
    //
    //
    void DefineHamiltonian3(int p, Complex coefficient);
    
    //
    //
    void DefineHamiltonian4(int p1, int p2, Complex coefficient);
    
    //
    //
    void DefineTensorHamiltonian();
};



#endif // DMRGCC_DMRG_COMPLEX_SPINFUL_BOSE_SQUARE_HAMILTONIAN_H_

