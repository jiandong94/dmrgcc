#include "dmrg/real_tensor_space.h"


RealTensorSpace::~RealTensorSpace()
{
    for(int i=0;i<num_site_pp_;++i)
    {
        for(int j=0;j<num_table_[i];++j)
        {
            delete[] quantum_table_[i][j];
        }
        delete[] quantum_table_[i];
    }

    delete[] quantum_table_;
    delete[] num_table_;

    for(int i=0;i<num_site_;++i) delete tensor_lattice_[i];

    delete[] tensor_lattice_;
}

RealTensorLattice* RealTensorSpace::get_tensor_lattice(int site)
{
    return tensor_lattice_[site];
}

void RealTensorSpace::DefineTensorSpaceCache(bool disk_cache, char *cache_name)
{
    disk_cache_ = disk_cache;
    strcpy(cache_name_, cache_name);
}

void RealTensorSpace::DefineTensorSpaceParameter(int max_block, int max_dim, 
        double canonical_precision, double noise_factor)
{
    max_block_ = max_block;
    max_dim_ = max_dim;
    canonical_precision_ = canonical_precision;
    noise_factor_ = noise_factor;
}

void RealTensorSpace::InitializeTensorSpace(int initial_block, int initial_dim)
{

}

void RealTensorSpace::PrintTensorSpace()
{
    cout << "===============TENSOR SPACE================" << endl;
    cout << "number of sites = " << num_site_ << endl;
    cout << "Quantum Table: " << endl;
    if(quantum_table_ == nullptr)
    {
        cout << "quantum table is nullptr!" << endl;
    }
    else
    {
        for(int i=0;i<num_site_pp_;++i)
        {
            cout << "----------Position="<< i << "----------" << endl;
            for(int k=0;k<num_quantum_;++k)
            {
                cout << "QN" << k << ": ";
                for(int j=0;j<num_table_[i];++j)
                {
                    cout << quantum_table_[i][j][k] << " ";
                }
                cout << endl;
            }
        }
    }
    cout << "Tensor Lattice: " << endl;
    for(int i=0;i<num_site;++i)
    {
        cout << "----------Site="<< i << "----------" << endl;
        tensor_lattice_[i]->PrintTensorLattice();
    }
}

void RealTensorSpace::WriteTensorSpace(const char* tensor_space_name)
{
    ofstream tensor_space_file;

    tensor_space_file.open(tensor_space_name, ios::binary|ios::out);

    WriteTensorSpace(tensor_space_file);

    tensor_space_file.close();
}

void RealTensorSpace::WriteTensorSpace(ofstream &tensor_space_file)
{
    tensor_space_file.write((char*) &num_site_, sizeof(int));
    tensor_space_file.write((char*) &num_site_pp_, sizeof(int));
    tensor_space_file.write((char*) &num_site_mm_, sizeof(int));

    tensor_space_file.write((char*) &num_quant_, sizeof(int));
    for(int i=0;i<num_site_pp_;++i)
        tensor_space_file.write((char*) &num_table_[i], sizeof(int));
    for(int i=0;i<num_site_pp_;++i)
        for(int j=0;j<num_table_[i];++j)
            for(int k=0;k<num_quantum_;++k)
                tensor_space_file.write((char*) &quantum_table_[i][j][k], sizeof(int));
    
    for(int i=0;i<num_site_;++i)
    {
        if(disk_cache_ == true) ResumeTensorSpace(i);
        tensor_lattice_[i]->WriteTensorLattice(tensor_space_file);
        if(disk_cache_ == true) ResetTensorSpace(i);
    }

    tensor_space_file.write((char*) &max_block_, sizeof(int));
    tensor_space_file.write((char*) &max_dim_, sizeof(int));
    tensor_space_file.write((char*) &canonical_precision_, sizeof(double));
    tensor_space_file.write((char*) &noise_factor_, sizeof(double));
}


void RealTensorSpace::ReadTensorSpace(const char* tensor_space_name)
{
    ifstream tensor_space_file;

    tensor_space_file.open(tensor_space_name, ios::binary|ios::in);

    ReadTensorSpace(tensor_space_file);

    tensor_space_file.close();
}

void RealTensorSpace::ReadTensorSpace(ifstream &tensor_space_file)
{
    for(int i=0;i<num_site_pp_;++i)
    {
        for(int j=0;j<num_table_[i];++j)
        {
            delete[] quantum_table_[i][j];
        }
        delete[] quantum_table_[i];
    }

    delete[] quantum_table_;
    delete[] num_table_;

    for(int i=0;i<num_site_;++i) delete tensor_lattice_[i];

    delete[] tensor_lattice_;


    tensor_space_file.read((char*) &num_site_, sizeof(int));
    tensor_space_file.read((char*) &num_site_pp_, sizeof(int));
    tensor_space_file.read((char*) &num_site_mm_, sizeof(int));

    tensor_space_file.read((char*) &num_quant_, sizeof(int));
    num_table_ = new int[num_quant];
    for(int i=0;i<num_site_pp_;++i)
        tensor_space_file.read((char*) &num_table_[i], sizeof(int));
    
    quantum_table_ = new int** [num_quant];
    for(int i=0;i<num_site_pp_;++i)
    {
        quantum_table_[i] = new int* [num_table_[i]];
        for(int j=0;j<num_table_[i];++j)
        {
            quantum_table_[i][j] = new int[num_quantum_];
            for(int k=0;k<num_quantum_;++k)
                tensor_space_file.read((char*) &quantum_table_[i][j][k], sizeof(int));
        }
    }

    tensor_lattice_ = new RealTensorLattice* [num_site_];
    for(int i=0;i<num_site_;++i)
    {
        tensor_lattice_[i] = new RealTensorLattice();
        tensor_lattice_[i]->ReadTensorLattice(tensor_space_file);
        if(disk_cache_ == true) 
        {
            RecordTensorSpace(i);
            ResetTensorSpace(i);
        }
    }

    tensor_space_file.read((char*) &max_block_, sizeof(int));
    tensor_space_file.read((char*) &max_dim_, sizeof(int));
    tensor_space_file.read((char*) &canonical_precision_, sizeof(double));
    tensor_space_file.read((char*) &noise_factor_, sizeof(double));
}

void RealTensorSpace::RecordTensorSpace(int site)
{
    char tensor_lattice_name[512];

    sprintf(tensor_lattice_name, "%s/SpaceDiskCacheFile/Tensor_Space_%d.dat", cache_name_, site);

    tensor_lattice_[site]->WriteTensorLattice(tensor_lattice_name);
}


void RealTensorSpace::ResumeTensorSpace(int site)
{
    char tensor_lattice_name[512];

    sprintf(tensor_lattice_name, "%s/SpaceDiskCacheFile/Tensor_Space_%d.dat", cache_name_, site);

    tensor_lattice_[site]->ReadTensorLattice(tensor_lattice_name);
}

void RealTensorSpace::ResetTensorSpace(int site)
{
    tensor_lattice_[site]->ResetTensorLattice();
}











