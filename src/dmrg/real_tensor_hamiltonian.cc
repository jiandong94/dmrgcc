#include "dmrg/real_tensor_hamiltonian.h"

RealTensorHamiltonian::~RealTensorHamiltonian()
{
    for(int i=0;i<num_site_pp_;++i)
    {
        for(int j=0;j<num_table_[i];++j) delete[] quantum_table_[i][j];
        delete[] quantum_table_[i];
    }
    delete[] quantum_table_;
    delete[] num_table_;

    for(int i=0;i<num_operator_;++i)
        delete basic_operator_[i];
    delete[] basic_operator_;
    
    for(int i=0;i<num_site_;++i)
        delete tensor_hamiltonian_[i];
    delete[] tensor_hamiltonian_;
}

RealTensorOperator* RealTensorHamiltonian::get_tensor_hamiltonian(int site)
{
    return tensor_hamiltonian_[site];
}

void RealTensorHamiltonian::PrintTensorHamiltonian()
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

    tensor_hamiltonian_file.open(tensor_hamiltonian_name, ios::binary|ios::out);

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
        tensor_hamiltonian_file.write((char*) &num_table_[i], sizeof(int));
    for(int i=0;i<num_site_pp_;++i)
        for(int j=0;j<num_table_[i];++j)
            for(int k=0;k<num_quantum_;++k)
                tensor_hamiltonian_file.write((char*) &quantum_table_[i][j][k], sizeof(int));

    for(int i=0;i<num_site_;++i)
        tensor_hamiltonian_[i]->WriteTensorOperator(tensor_hamiltonian_file);
}

void RealTensorHamiltonian::ReadTensorHamiltonian(const char* tensor_hamiltonian_name)
{
    ifstream tensor_hamiltonian_file;

    tensor_hamiltonian_file.open(tensor_hamiltonian_name, ios::binary|ios::in);

    ReadTensorHamiltonian(tensor_hamiltonian_file);

    tensor_hamiltonian_file.close();
}

void RealTensorHamiltonian::ReadTensorHamiltonian(ifstream &tensor_hamiltonian_file)
{
    for(int i=0;i<num_site_pp_;++i)
    {
        for(int j=0;j<num_table_[i];++j) delete[] quantum_table_[i][j];
        delete[] quantum_table_[i];
    }
    delete[] quantum_table_;
    delete[] num_table_;

    for(int i=0;i<num_operator_;++i)
        delete basic_operator_[i];
    delete[] basic_operator_;
    
    for(int i=0;i<num_site_;++i)
        delete tensor_hamiltonian_[i];
    delete[] tensor_hamiltonian_;

    tensor_hamiltonian_file.read((char*) &num_site_, sizeof(int));
    tensor_hamiltonian_file.read((char*) &num_site_pp_, sizeof(int));
    tensor_hamiltonian_file.read((char*) &num_site_mm_, sizeof(int));
    
    tensor_hamiltonian_file.read((char*) &num_operator_, sizeof(int));
    basic_operator_ = new RealMatrix*[num_operator_];
    for(int i=0;i<num_operator_;++i)
        basic_operator_[i]->ReadMatrix(tensor_hamiltonian_file);

    tensor_hamiltonian_file.read((char*) &num_quantum_, sizeof(int));
    num_table_ = new int[num_site_pp_];
    for(int i=0;i<num_site_pp_;++i)
        tensor_hamiltonian_file.read((char*) &num_table_[i], sizeof(int));
    quantum_table_ = new int** [num_site_pp_];
    for(int i=0;i<num_site_pp_;++i)
    {
        quantum_table_[i] = new int* [num_table_[i]];
        for(int j=0;j<num_table_[i];++j)
        {
            quantum_table_[i][j] = new int[num_quantum_];
            for(int k=0;k<num_quantum_;++k)
                tensor_hamiltonian_file.read((char*) &quantum_table_[i][j][k], sizeof(int));
        }
    }
}

void RealTensorHamiltonian::ExpanTensorHamiltonian(int site, int expan_operator_index, int expan_coefficient, 
        int* expan_table)
{
    int expan_num_table, **expan_quantum_table;
    int leigh = -1;
    if(site == 0) leigh = 0;
    if(site == num_site_mm_) leigh = 1;

    // expan tensor operator
    tensor_hamiltonian_[site]->ExpanTensorOperator(basic_operator_, leigh, 
            expan_operator_index, expan_coefficient);

    // expan quantum table
    if(site < num_site_mm_)
    {
        expan_num_table = num_table_[site+1]+1;
        expan_quantum_table = new int* [expan_num_table];
        for(int i=0;i<expan_num_table;++i)
            expan_quantum_table[i] = new int[num_quantum_];
        for(int j=0;j<num_quantum_;++j)
        {
            for(int i=0;i<num_table_[site+1];++i)
            { 
                // pervious quantum table
                expan_quantum_table[i][j] = quantum_table_[site+1][i][j];
            }
            // previous quantum table + expaned quantum_table
            expan_quantum_table[num_table_[site+1]][j] = expan_table[j];
        }

        for(int i=0;i<num_table_[site+1];++i) delete[] quantum_table_[site+1][i];
        delete[] quantum_table_[site+1];

        num_table_[site+1] = expan_num_table;
        quantum_table_[site+1] = expan_quantum_table;
    }
}

void ParallelTensorHamiltonian()
{
    RealMatrix* transfer_tensor;
    int num_unparallel, *position_unparallel;

    for(int i=0;i<num_site_mm_;++i)
    {
        tensor_hamiltonian_[i]->LeftParallelTensorOperator(num_unparallel, position_unparallel, 
                                                           transfer_tensor);
        tensor_hamiltonian_[i+1]->LeftMergeTensorOperator(num_unparallel, transfer_tensor);
        ParallelQuantumTable(i+1, num_unparallel, position_unparallel);

        delete[] position_unparallel;
        delete transfer_tensor;
    }
    for(int i=num_site_mm_;i>0;--i)
    {
        tensor_hamiltonian_[i]->RightParallelTensorOperator(num_unparallel, position_unparallel, 
                                                           transfer_tensor);
        tensor_hamiltonian_[i-1]->RightMergeTensorOperator(num_unparallel, transfer_tensor);
        ParallelQuantumTable(i, num_unparallel, position_unparallel);

        delete[] position_unparallel;
        delete transfer_tensor;
    }
}

void RealTensorHamiltonian::DefineQuantumTable(int num_quantum)
{
    num_quantum_ = num_quantum;
    num_table_ = new int[num_site_pp_];
    quantum_table_ = new int** [num_site_pp_];
    for(int i=0;i<num_site_pp_;++i)
    {
        num_table_[i] = 0;
        quantum_table_[i] = nullptr;
    }
}

void RealTensorHamiltonian::ParallelQuantumTable(int site, int num_unparallel, 
        int* position_unparallel)
{
    int **result_quantum_table;

    result_quantum_table = new int* [num_unparallel];
    for(int i=0;i<num_unparallel;++i)
    {
        result_quantum_table = new int[num_quantum_];
        for(int j=0;j<num_quantum_;++j)
        {
            result_quantum_table[i][j] = quantum_table_[site][position_unparallel[i]][j]
        }
    }

    for(int i=0;i<num_table_[site];++i)
        delete[] quantum_table_[site][i];
    delete[] quantum_table_[site];

    num_table_[site] = num_unparallel;
    quantum_table_[site] = result_quantum_table;
}
