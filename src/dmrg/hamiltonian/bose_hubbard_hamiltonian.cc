#include "dmrg/hamiltonian/bose_hubbard_hamiltonian.h"


BoseHubbardHamiltonian::BoseHubbardHamiltonian(int num_site_x, int num_site_y, int physics_dim, 
        double flux_value, double hop_x, double hop_y, double inter_value, 
        double chemical_value=0, bool period_x = false, bool period_y = false)
{
    num_site_x_ = num_site_x;
    num_site_y_ = num_site_y;
    physics_dim_ = physics_dim;

    flux_value_ = flux_value;
    hop_x_ = hop_x;
    hop_y_ = hop_y;
    inter_value_ = inter_value;
    chemical_value_ = chemical_value;
    
    period_x_ = period_x;
    period_y_ = period_y;


    num_site_ = num_site_x*num_site_y;
    num_site_pp_ = num_site_+1;
    num_site_mm_ = num_site_-1;

    DefineBasicTensor();
    DefineQuantumTable();

    tensor_hamiltonian_ = new RealTensorOperator* [num_site_];
    for(int i=0;i<num_site;++i)
    {
        tensor_hamiltonian_ = new RealTensorOperator(physics_dim_);
    }
}

BoseHubbardHamiltonian::~BoseHubbardHamiltonian()
{

}

BoseHubbarHamiltonian::DefineBasicTensor()
{
    // 0=I, 1=a^+, 2=a^-, 3=n, 4={\mu}n+{U}n(n-1)
    num_operator_ = 5;

    basic_tensor_ = new RealMatrix* [num_operator_];
    for(int i=0;i<num_operator_;++i)
    {
        basic_tensor_ = new RealMatrix(physics_dim_, physics_dim_);
    }

    for(int i=0;i<physics_dim_;++i)
    {
        basic_tensor_[0]->set_matrix_element(i, i, 1.0);
        basic_tensor_[3]->set_matrix_element(i, i, i);
        basic_tensor_[4]->set_matrix_element(i, i, chemical_value_*i+inter_value_*i*(i-1));
    }

    for(int i=1;i<physics_dim_;++i)
    {
        basic_tensor_[1]->set_matrix_element(i, i-1, sqrt(i));
        basic_tensor_[2]->set_matrix_element(i-1, i, sqrt(i));      
    }    
}

