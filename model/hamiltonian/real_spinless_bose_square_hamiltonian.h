#ifndef DMRGCC_DMRG_REAL_SPINLESS_BOSE_SQUARE_HAMILTONIAN_H_
#define DMRGCC_DMRG_REAL_SPINLESS_BOSE_SQUARE_HAMILTONIAN_H_

#include "dmrg/realdmrg/real_tensor_hamiltonian.h"
class RealSpinlessBoseSquareHamiltonian : public RealTensorHamiltonian
{
    protected:
    
    int num_site_x_;
    int num_site_y_;
    int physics_dim_;

    double flux_value_;
    double hop_x_;
    double hop_y_;
    double inter_value_;
    double chemical_value_;
    
    bool period_x_;
    bool period_y_;

    public:

    //
    //
    RealSpinlessBoseSquareHamiltonian(int num_site_x, int num_site_y, int physics_dim, double flux_value, 
            double hop_x, double hop_y, double inter_value, double chemical_value=0, 
            bool period_x=false, bool period_y=false);

    //
    //
    ~RealSpinlessBoseSquareHamiltonian();

    //
    //
    void DefineBasicTensor();

    //
    //
    void DefineHamiltonian1(int p1, int p2, double coefficient);

    //
    //
    void DefineHamiltonian2(int p1, int p2, double coefficient);
    
    //
    //
    void DefineHamiltonian3(int p, double coefficient);
    
    //
    //
    void DefineTensorHamiltonian();
};



#endif // DMRGCC_DMRG_REAL_SPINLESS_BOSE_SQUARE_HAMILTONIAN_H_

