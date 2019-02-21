#ifndef DMRGCC_DMRG_REAL_TENSOR_HAMILTONIAN_H_
#define DMRGCC_DMRG_REAL_TENSOR_HAMILTONIAN_H_

#include "dmrg/realdmrg/real_tensor_operator.h"

class RealTensorHamiltonian
{
    friend class RealTensorContraction;
    friend class RealTensorNetwork;
    protected:

    int num_site_;
    int num_site_pp_; // number of sites++
    int num_site_mm_; // number of sites--

    int num_operator_;

    RealMatrix** basic_operator_;

    int num_quantum_;
    int *num_table_;
    int ***quantum_table_;

    RealTensorOperator **tensor_hamiltonian_;

    public:

    //
    //
    virtual ~RealTensorHamiltonian();

    //
    //
    RealTensorOperator* get_tensor_hamiltonian(int site);

    //
    //
    void PrintTensorHamiltonian();

    //
    //
    void WriteTensorHamiltonian(const char* tensor_hamiltonian_name);

    //
    //
    void WriteTensorHamiltonian(ofstream &tensor_hamiltonian_file);

    //
    //
    void ReadTensorHamiltonian(const char* tensor_hamiltonian_name);

    //
    //
    void ReadTensorHamiltonian(ifstream &tensor_hamiltonian_file);

    //
    //
    void ExpanTensorHamiltonian(int site, int expan_operator_index, double expan_coefficient, 
            int* expan_table);

    //
    //
    void ParallelTensorHamiltonian();

    //
    //
    virtual void DefineTensorHamiltonian();

    protected:

    //
    //
    void DefineQuantumTable(int num_quantum);

    //
    //
    void ParallelQuantumTable(int site, int num_unparallel, int* position_unparallel);
};



#endif // DMRGCC_DMRG_REAL_TENSOR_HAMILTONIAN_H_
