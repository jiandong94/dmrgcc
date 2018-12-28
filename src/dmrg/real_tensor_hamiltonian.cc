#include "dmrg/real_tensor_hamiltonian.h"

RealTensorHamiltonian::~RealTensorHamiltonian()
{
    for(int i=0;i<num_site_pp_;++i)
    {
        for(int j=0;j<num_table_[i];++j) delete[] quantum_table_[i][j];
        delete[] quantum_table_[i]
    }
    delete[] quantum_table_;
    delete[] num_table_;

    for(int i=0;i<num_operator_;++i)
        delete basic_operator_[i];
    delete[] basic_operator;
    
    for(int i=0;i<num_site_;++i)
        delete tensor_hamiltonian_[i];
    delete[] tensor_hamiltonian_;
}

RealTensorOperator* RealTensorHamiltonian::get_tensor_hamiltonian(site)
{
    return tensor_hamiltonian_[site];
}

void RealTensorHamiltonian::PrintHamiltonian()
{
    cout << "==============================" << endl;
    cout << "TensorHamiltonian" << endl;
    for(int i=0;i<num_site_;++i)
    {
        cout << "----------site = " << i << "----------" << endl;
        tensor_hamiltonian_[i]->PrintTensorOperator();
    }
}

void RealTensorHamiltonian::WriteTensorHamiltonian(const char* tensor_hamiltonian_name)
{
    ofstream tensor_hamiltonian_file;

    tenosr_hamiltonian_file.open(tensor_hamiltonian_name, ios::binary|ios::out);

    WriteTensorHamiltonian(tensor_hamiltonian_file);

    tensor_hamiltonian_file.close();
}


void RealTensorHamiltonian::WriteTensorHamiltonian(ofstream &tensor_hamiltonian_file)
{
    tensor_hamiltonian_file.write((char*) &num_site_, sizeof(int));
    tensor_hamiltonian_file.write((char*) &num_site_pp_, sizeof(int));
    tensor_hamiltonian_file.write((char*) &num_site_mm_, sizeof(int));
    
    tensor_hamiltonian_file.write((char*) &num_operator_, sizeof(int));
    for(int i=0;i<num_operator_;++i)
        basic_operator_[i]->WriteMatrix(tensor_hamiltonian_file);

    tensor_hamiltonian_file.write((char*) &num_quantum_, sizeof(int));
    for(int i=0;i<num_site_pp_;++i)
        tensor_hamiltonian_file.write((char*) &num_table_[i]_, sizeof(int));
    for(int i=0;i<num_site_pp_;++i)
        for(int j=0;j<num_table_[i];++j)
            for(int k=0;k<num_quantum_;++k)
                tensor_hamiltonian_file.write((char*) &quantum_table_[i][j][k]_, sizeof(int));

    for(int i=0;i<num_site_;++i)
        tensor_hamiltonian_[i]->WriteTensorOperator[tensor_hamiltonian_file];
}

void RealTensorHamiltonian::ReadTensorHamiltonian(const char* tensor_hamiltonian_name)
{
    ifstream tensor_hamiltonian_file;

    tenosr_hamiltonian_file.open(tensor_hamiltonian_name, ios::binary|ios::in);

    ReadTensorHamiltonian(tensor_hamiltonian_file);

    tensor_hamiltonian_file.close();
}

void RealTensorHamiltonian::ReadTensorHamiltonian(ifstream &tensor_hamiltonian_file)
{
    for(int i=0;i<num_site_pp_;++i)
    {
        for(int j=0;j<num_table_[i];++j) delete[] quantum_table_[i][j];
        delete[] quantum_table_[i]
    }
    delete[] quantum_table_;
    delete[] num_table_;

    for(int i=0;i<num_operator_;++i)
        delete basic_operator_[i];
    delete[] basic_operator;
    
    for(int i=0;i<num_site_;++i)
        delete tensor_hamiltonian_[i];
    delete[] tensor_hamiltonian_;



}
