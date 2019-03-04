#include "model/hamiltonian/real_spinless_bose_square_hamiltonian.h"


RealSpinlessBoseSquareHamiltonian::RealSpinlessBoseSquareHamiltonian(int num_site_x, int num_site_y, int physics_dim, 
        double flux_value, double hop_x, double hop_y, double inter_value, 
        double chemical_value, bool period_x, bool period_y)
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
    DefineQuantumTable(1);

    tensor_hamiltonian_ = new RealTensorOperator* [num_site_];
    for(int i=0;i<num_site_;++i)
    {
        tensor_hamiltonian_[i] = new RealTensorOperator(physics_dim_);
    }
}

RealSpinlessBoseSquareHamiltonian::~RealSpinlessBoseSquareHamiltonian()
{

}

void RealSpinlessBoseSquareHamiltonian::DefineBasicTensor()
{
    // 0=I, 1=a^+, 2=a^-, 3=n, 4={\mu}n+{U}n(n-1)
    num_operator_ = 5;

    basic_operator_ = new RealMatrix* [num_operator_];
    for(int i=0;i<num_operator_;++i)
    {
        basic_operator_[i] = new RealMatrix(physics_dim_, physics_dim_);
    }

    for(int i=0;i<physics_dim_;++i)
    {
        basic_operator_[0]->set_matrix_element(i, i, 1.0);
        basic_operator_[3]->set_matrix_element(i, i, i);
        basic_operator_[4]->set_matrix_element(i, i, chemical_value_*i+inter_value_*i*(i-1));
    }

    for(int i=1;i<physics_dim_;++i)
    {
        basic_operator_[1]->set_matrix_element(i, i-1, sqrt(i));
        basic_operator_[2]->set_matrix_element(i-1, i, sqrt(i));      
    }    
}

// a^+ a^-
void RealSpinlessBoseSquareHamiltonian::DefineHamiltonian1(int p1, int p2, double coefficient)
{
    double operator_coefficient;
    int operator_index, operator_table[1];

    if(p1>=p2)
    {
        error("position1 >= position2 in DefineHamiltonian1");
    }

    for(int i=0;i<num_site_;++i)
    {
        operator_index = 0; // I
        operator_coefficient = 1.0;
        operator_table[0] = 0;

        if(i == p1)
        {
            operator_index = 1; // a^+
            operator_coefficient = 1.0;
            operator_table[0] = 1;
        }

        if(i>p1 && i<p2)
        {
            operator_index = 0; // I
            operator_coefficient = 1.0;
            operator_table[0] = 1;
        }
        if(i == p2)
        {
            operator_index = 2; // a^-
            operator_coefficient = coefficient;
            operator_table[0] = 0;
        }
        
        ExpanTensorHamiltonian(i, operator_index, operator_coefficient, operator_table);
    }
}

// a^- a^+
void RealSpinlessBoseSquareHamiltonian::DefineHamiltonian2(int p1, int p2, double coefficient)
{
    double operator_coefficient;
    int operator_index, operator_table[1];

    if(p1>=p2)
    {
        error("position1 >= position2 in DefineHamiltonian2");
    }

    for(int i=0;i<num_site_;++i)
    {
        operator_index = 0; // I
        operator_coefficient = 1.0;
        operator_table[0] = 0;

        if(i == p1)
        {
            operator_index = 2; // a^-
            operator_coefficient = 1.0;
            operator_table[0] = -1;
        }

        if(i>p1 && i<p2)
        {
            operator_index = 0; // I
            operator_coefficient = 1.0;
            operator_table[0] = -1;
        }
        if(i == p2)
        {
            operator_index = 1; // a^+
            operator_coefficient = coefficient;
            operator_table[0] = 0;
        }
        
        ExpanTensorHamiltonian(i, operator_index, operator_coefficient, operator_table);
    }
}

// {\mu}n+{U}n(n-1)
void RealSpinlessBoseSquareHamiltonian::DefineHamiltonian3(int p, double coefficient)
{
    double operator_coefficient;
    int operator_index, operator_table[1];

    for(int i=0;i<num_site_;++i)
    {
        operator_index = 0; // I
        operator_coefficient = 1.0;
        operator_table[0] = 0;

        if(i == p)
        {
            operator_index = 4; // a^+
            operator_coefficient = coefficient;
            operator_table[0] = 0;
        }

        ExpanTensorHamiltonian(i, operator_index, operator_coefficient, operator_table);
    }
}

void RealSpinlessBoseSquareHamiltonian::DefineTensorHamiltonian()
{
    int p[2];

    // quantum table
    num_table_[0] = 1;
    num_table_[num_site_] = 1;

    quantum_table_[0] = new int* [1];
    quantum_table_[0][0] = new int[1];
    quantum_table_[0][0][0] = 0;

    quantum_table_[num_site_] = new int* [1];
    quantum_table_[num_site_][0] = new int[1];
    quantum_table_[num_site_][0][0] = 0;

    for(int i=0;i<num_site_x_;++i) for(int j=0;j<num_site_y_;++j)
    {
        p[0] = j+i*num_site_y_;
        
        if(hop_x_ != 0)
        {
            p[1] = j+((i+1)%num_site_x_)*num_site_y_;
            if(i+1 < num_site_x_)
            {
                DefineHamiltonian1(p[0], p[1], hop_x_);
                DefineHamiltonian2(p[0], p[1], hop_x_);
            }
            else if(num_site_x_>2 && period_x_==true)
            {
                DefineHamiltonian1(p[1], p[0], hop_x_);
                DefineHamiltonian2(p[1], p[0], hop_x_);
            }
        }
        if(hop_y_ != 0)
        {
            p[1] = (j+1)%num_site_y_+i*num_site_y_;
            if(j+1 < num_site_y_)
            {
                DefineHamiltonian1(p[0], p[1], hop_y_);
                DefineHamiltonian2(p[0], p[1], hop_y_);
            }
            else if(num_site_y_>2 && period_y_==true)
            {
                DefineHamiltonian1(p[1], p[0], hop_y_);
                DefineHamiltonian2(p[1], p[0], hop_y_);
            }
        }

        if(chemical_value_!=0 || inter_value_!=0)
        {
            DefineHamiltonian3(p[0], 1.0);
        }
    
        ParallelTensorHamiltonian();
    }
    //PrintTensorHamiltonian();
}
