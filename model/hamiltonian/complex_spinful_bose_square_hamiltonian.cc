#include "model/hamiltonian/complex_spinful_bose_square_hamiltonian.h"


ComplexSpinfulBoseSquareHamiltonian::ComplexSpinfulBoseSquareHamiltonian(int num_site_x, int num_site_y, int physics_dim, 
        double* flux, double* hop_x, double* hop_y, double* inter_value, 
        double* chemical_value, bool period_x, bool period_y)
{
    num_site_x_ = num_site_x;
    num_site_y_ = num_site_y;
    physics_dim_ = physics_dim;

    for(int t=0;t<2;++t)
    {
        flux_[t] = flux[t];
        hop_x_[t] = hop_x[t];
        hop_y_[t] = hop_y[t];
        chemical_value_[t] = chemical_value[t];
    }
    for(int t=0;t<3;++t)
        inter_value_[t] = inter_value[t];
    
    period_x_ = period_x;
    period_y_ = period_y;


    num_site_ = 2*num_site_x*num_site_y;
    num_site_pp_ = num_site_+1;
    num_site_mm_ = num_site_-1;
    
    DefineBasicTensor();
    DefineQuantumTable(2);

    tensor_hamiltonian_ = new ComplexTensorOperator* [num_site_];
    for(int i=0;i<num_site_;++i)
    {
        tensor_hamiltonian_[i] = new ComplexTensorOperator(physics_dim_);
    }
}

ComplexSpinfulBoseSquareHamiltonian::~ComplexSpinfulBoseSquareHamiltonian()
{

}

void ComplexSpinfulBoseSquareHamiltonian::DefineBasicTensor()
{
    // 0=I, 1=a^+, 2=a^-, 3=n, 4={\mu_0}n+0.5*{U_0}n(n-1), 5={\mu_1}n+0.5*{U_1}n(n-1)
    num_operator_ = 6;

    basic_operator_ = new ComplexMatrix* [num_operator_];
    for(int i=0;i<num_operator_;++i)
    {
        basic_operator_[i] = new ComplexMatrix(physics_dim_, physics_dim_);
    }

    for(int i=0;i<physics_dim_;++i)
    {
        basic_operator_[0]->set_matrix_element(i, i, 1.0);
        basic_operator_[3]->set_matrix_element(i, i, i);
        basic_operator_[4]->set_matrix_element(i, i, chemical_value_[0]*i+0.5*inter_value_[0]*i*(i-1));
        basic_operator_[5]->set_matrix_element(i, i, chemical_value_[0]*i+0.5*inter_value_[1]*i*(i-1));
    }

    for(int i=1;i<physics_dim_;++i)
    {
        basic_operator_[1]->set_matrix_element(i, i-1, sqrt(i));
        basic_operator_[2]->set_matrix_element(i-1, i, sqrt(i));      
    }    
}

// a^+ a^-
void ComplexSpinfulBoseSquareHamiltonian::DefineHamiltonian1(int p1, int p2, Complex coefficient)
{
    Complex operator_coefficient;
    int operator_index, operator_table[2], qn_operator[2];

    if(p1>=p2)
    {
        error("position1 >= position2 in DefineHamiltonian1");
    }

    if(p1%2 == 0)
    {
        qn_operator[0] = 1;
        qn_operator[1] = 0;
    }
    else
    {
        qn_operator[0] = 0;
        qn_operator[1] = 1;
    }

    for(int i=0;i<num_site_;++i)
    {
        operator_index = 0; // I
        operator_coefficient = 1.0;
        operator_table[0] = 0;
        operator_table[1] = 0;

        if(i == p1)
        {
            operator_index = 1; // a^+
            operator_coefficient = 1.0;
            operator_table[0] = qn_operator[0];
            operator_table[1] = qn_operator[1];
        }

        if(i>p1 && i<p2)
        {
            operator_index = 0; // I
            operator_coefficient = 1.0;
            operator_table[0] = qn_operator[0];
            operator_table[1] = qn_operator[1];
        }
        if(i == p2)
        {
            operator_index = 2; // a^-
            operator_coefficient = coefficient;
            operator_table[0] = 0;
            operator_table[1] = 0;
        }
        
        ExpanTensorHamiltonian(i, operator_index, operator_coefficient, operator_table);
    }
}

// a^- a^+
void ComplexSpinfulBoseSquareHamiltonian::DefineHamiltonian2(int p1, int p2, Complex coefficient)
{
    Complex operator_coefficient;
    int operator_index, operator_table[2], qn_operator[2];

    if(p1>=p2)
    {
        error("position1 >= position2 in DefineHamiltonian2");
    }
    
    if(p1%2 == 0)
    {
        qn_operator[0] = -1;
        qn_operator[1] = 0;
    }
    else
    {
        qn_operator[0] = 0;
        qn_operator[1] = -1;
    }

    for(int i=0;i<num_site_;++i)
    {
        operator_index = 0; // I
        operator_coefficient = 1.0;
        operator_table[0] = 0;
        operator_table[1] = 0;

        if(i == p1)
        {
            operator_index = 2; // a^-
            operator_coefficient = 1.0;
            operator_table[0] = qn_operator[0];
            operator_table[1] = qn_operator[1];
        }

        if(i>p1 && i<p2)
        {
            operator_index = 0; // I
            operator_coefficient = 1.0;
            operator_table[0] = qn_operator[0];
            operator_table[1] = qn_operator[1];
        }
        if(i == p2)
        {
            operator_index = 1; // a^+
            operator_coefficient = coefficient;
            operator_table[0] = 0;
            operator_table[1] = 0;
        }
        
        ExpanTensorHamiltonian(i, operator_index, operator_coefficient, operator_table);
    }
}

// {\mu}n+0.5*{U}n(n-1)
void ComplexSpinfulBoseSquareHamiltonian::DefineHamiltonian3(int p, Complex coefficient)
{
    Complex operator_coefficient;
    int operator_index, operator_table[2];

    for(int i=0;i<num_site_;++i)
    {
        operator_index = 0; // I
        operator_coefficient = 1.0;
        operator_table[0] = 0;
        operator_table[1] = 0;

        if(i == p)
        {
            if(p%2 == 0)
                operator_index = 4; // a^+
            else
                operator_index = 5;
            operator_coefficient = coefficient;
            operator_table[0] = 0;
            operator_table[1] = 0;
        }

        ExpanTensorHamiltonian(i, operator_index, operator_coefficient, operator_table);
    }
}

// 
void ComplexSpinfulBoseSquareHamiltonian::DefineHamiltonian4(int p1, int p2, Complex coefficient)
{
    Complex operator_coefficient;
    int operator_index, operator_table[2];
    for(int i=0;i<num_site_;++i)
    {
        operator_index = 0; // I
        operator_coefficient = 1.0;
        operator_table[0] = 0;
        operator_table[1] = 0;

        if(i == p1)
        {
            operator_index = 3; // n
            operator_coefficient = 1.0;
            operator_table[0] = 0;
            operator_table[1] = 0;
        }
        
        if(i == p2)
        {
            operator_index = 3; // n
            operator_coefficient = coefficient;
            operator_table[0] = 0;
            operator_table[1] = 0;
        }
        ExpanTensorHamiltonian(i, operator_index, operator_coefficient, operator_table);
    }
}

void ComplexSpinfulBoseSquareHamiltonian::DefineTensorHamiltonian()
{
    int p[2];
    Complex phase;

    // quantum table
    num_table_[0] = 1;
    num_table_[num_site_] = 1;

    // first site
    quantum_table_[0] = new int* [1];
    quantum_table_[0][0] = new int[2];
    quantum_table_[0][0][0] = 0;
    quantum_table_[0][0][1] = 0;

    // last site
    quantum_table_[num_site_] = new int* [1];
    quantum_table_[num_site_][0] = new int[2];
    quantum_table_[num_site_][0][0] = 0;
    quantum_table_[num_site_][0][1] = 0;

    for(int i=0;i<num_site_x_;++i) for(int j=0;j<num_site_y_;++j)
    {
        p[0] = j+i*num_site_y_;
        
        if(hop_x_[0] != 0.0)
        {
            p[1] = j+((i+1)%num_site_x_)*num_site_y_;
            if(i+1 < num_site_x_)
            {
                DefineHamiltonian1(2*p[0], 2*p[1], hop_x_[0]);
                DefineHamiltonian2(2*p[0], 2*p[1], hop_x_[0]);
            }
            else if(num_site_x_>2 && period_x_==true)
            {
                DefineHamiltonian1(2*p[1], 2*p[0], hop_x_[0]);
                DefineHamiltonian2(2*p[1], 2*p[0], hop_x_[0]);
            }
        }

        if(hop_x_[1] != 0.0)
        {
            p[1] = j+((i+1)%num_site_x_)*num_site_y_;
            if(i+1 < num_site_x_)
            {
                DefineHamiltonian1(2*p[0]+1, 2*p[1]+1, hop_x_[1]);
                DefineHamiltonian2(2*p[0]+1, 2*p[1]+1, hop_x_[1]);
            }
            else if(num_site_x_>2 && period_x_==true)
            {
                DefineHamiltonian1(2*p[1]+1, 2*p[0]+1, hop_x_[1]);
                DefineHamiltonian2(2*p[1]+1, 2*p[0]+1, hop_x_[1]);
            }
        }
        
        if(hop_y_[0] != 0.0)
        {
            p[1] = (j+1)%num_site_y_+i*num_site_y_;

            phase = Phase(2*M_PI*(i)*flux_[0]);

            if(j+1 < num_site_y_)
            {
                DefineHamiltonian1(2*p[0], 2*p[1], Conj(phase)*hop_y_[0]);
                DefineHamiltonian2(2*p[0], 2*p[1], phase*hop_y_[0]);
            }
            else if(num_site_y_>2 && period_y_==true)
            {
                DefineHamiltonian2(2*p[1], 2*p[0], Conj(phase)*hop_y_[0]);
                DefineHamiltonian1(2*p[1], 2*p[0], phase*hop_y_[0]);
            }
        }

        if(hop_y_[1] != 0.0)
        {
            p[1] = (j+1)%num_site_y_+i*num_site_y_;

            phase = Phase(2*M_PI*(i)*flux_[1]);

            if(j+1 < num_site_y_)
            {
                DefineHamiltonian1(2*p[0]+1, 2*p[1]+1, Conj(phase)*hop_y_[1]);
                DefineHamiltonian2(2*p[0]+1, 2*p[1]+1, phase*hop_y_[1]);
            }
            else if(num_site_y_>2 && period_y_==true)
            {
                DefineHamiltonian2(2*p[1]+1, 2*p[0]+1, Conj(phase)*hop_y_[1]);
                DefineHamiltonian1(2*p[1]+1, 2*p[0]+1, phase*hop_y_[1]);
            }
        }

        if(chemical_value_[0]!=0.0 || inter_value_[0]!=0.0)
        {
            DefineHamiltonian3(2*p[0], 1.0);
        }
    
        if(chemical_value_[1]!=0.0 || inter_value_[1]!=0.0)
        {
            DefineHamiltonian3(2*p[0]+1, 1.0);
        }
        
        if(inter_value_[2]!=0.0)
        {
            DefineHamiltonian4(2*p[0]+1, 2*p[0]+1, inter_value_[2]);
        }
        ParallelTensorHamiltonian();
    }
    //PrintTensorHamiltonian();
}
