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
    int num_left_block, num_right_block, *left_block, *right_block, *left_dim, *right_dim, 
        num_block, *left_index, *right_index, *physics_index, row, column;

    cout << "Initial Block = " << initial_block << ", Intial Dimension = " << initial_dim << endl;
    max_block_ = initial_block;
    max_dim_ = -1;
    noise_factor = 10;
    canonical_precision = 1e-30;

    // initial lattice block
    for(int i=0;i<num_site_;++i)
    {
        ComputeLatticeBlock(i, num_left_block, num_right_block, left_block, right_block);
        
        left_dim = new int[num_left_block];
        for(int l=0;l<num_left_block;++l)
        {
            if(i == 0) left_dim[l] = 1;
            else left_dim[l] = initial_dim;
        }
        
        right_dim = new int[num_right_block];
        for(int r=0;r<num_right_block;++r)
        {
            if(i == num_site_mm_) right_dim[r] = 1;
            else left_dim[r] = initial_dim;
        }

        tensor_lattice_[i]->DefineTensorLattice(num_left_block, num_right_block, left_block, 
                right_block, left_dim, right_dim);

        delete[] left_block;
        delete[] right_block;
        delete[] left_dim;
        delete[] right_dim;
    }

    // initial ket tensor
    // left canonical
    // disk cache
    for(int i=0;i<num_site_;++i)
    {
        ComputeTensorIndex(i, num_block, left_index, right_index, physics_index);

        row = initial_dim;
        column = initial_dim;
        
        if(i == 0) row = 1;
        if(i == num_site_mm_) column = 1;
        tensor_lattice_[i]->DefineTensorLattice(num_block, left_index, right_index, 
                physics_index, row, column);

        delete[] left_index;
        delete[] right_index;
        delete[] physics_index;

        CanonicalTensorSpace(0, i);
        if(i > 0) MergeTensorSpace(0, i);

        if(disk_cache_ = true)
        {
            RecordTensorSpace(i);
            ResetTensorSpace(i);
        }
    }
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

void RealTensorSpace::CanonicalTensorSpace(int leigh, int site)
{
    if(leigh==0) 
        tensor_lattice_[site]->LeftCanonicalTensorLattice(max_dim_, canonical_precision);
    else if(leigh==1)
        tensor_lattice_[site]->RightCanonicalTensorLattice(max_dim_, canonical_precision);
}

void RealTensorSpace::MergeTensorSpace(int leigh, int site)
{
    if(leigh==0 && site>0)
        tensor_lattice_[site]->LeftMergeTensorLattice(tensor_lattice_[site-1]);
    else if(leigh==1 && site<num_site_mm_)
        tensor_lattice_[site]->RightMergeTensorLattice(tensor_lattice_[site+1]);
}

void RealTensorSpace::ComputeLatticeBlock(int site, int &num_left_block, int &num_right_block, 
        int* &left_block, int* &right_block)
{
    int right_num_table, **left_quantum_table, **right_quantum_table, *tmp_left_block, 
        info;
    right_num_table = num_table_[site+1];
    left_quantum_table = quantum_table_[site];
    right_quantum_table = quantum_table_[site+1];
    
    // compute the left block from the previous site
    if(site == 0)
    {
        num_left_block = 0;
        left_block = new int[1];
        left_block[0] = 0;
    }
    else
    {
        num_left_block = tensor_lattice_[site]->get_num_left_block();
        left_block = new int[num_left_block];
        tmp_left_block = tensor_lattice_[site]->get_left_block();
        for(int i=0;i<num_left_block;++i)
            left_block[i] = tmp_left_block[i];
    }
    
    // compute number of right block
    num_right_block = 0;
    for(int r=0;r<right_num_table;++r)
    {
        for(int l=0;l<num_left_block;++l)
        {
            // left_block[l], not l
            info = CheckQuantumTable(site, left_quantum_table[left_block[l]], right_quantum_table[r]);
            if(info != -1)
            {
                num_right_block++;
                break;
            }
        }
        if(num_right_block == max_block)
            break;
    }

    // compute right block
    right_block = new int[num_right_block];
    int k = 0;
    for(int r=0;r<right_num_table;++r)
    {
        for(int l=0;l<num_left_block;++l)
        {
            // left_block[l], not l
            info = CheckQuantumTable(site, left_quantum_table[left_block[l]], right_quantum_table[r]);
            if(info != -1)
            {
                right_block[k] = r;
                k++;
                break;
            }
        }
        if(k == max_block)
            break;
    }
}

void RealTensorSpace::ComputeTensorIndex(int site, int &num_block, int* &left_index, int* &right_index, 
                         int* &physics_index)
{
    int num_left_block, num_right_block, *left_block, *right_block, **left_quantum_table, 
        **right_quantum_table, info;
    
    num_left_block = tensor_lattice_[site]->get_num_left_block();
    num_right_block = tensor_lattice_[site]->get_num_right_block();
    left_block = tensor_lattice_[site]->get_left_block();
    right_block = tensor_lattice_[site]->get_right_block();

    left_quantum_table = quantum_table_[site];
    right_quantum_table = quantum_table_[site+1];

    num_block = 0;
    for(int i=0;i<num_left_block;++i) for(int j=0;j<num_right_block;++j)
    {   
        info = CheckQuantumTable(site, left_quantum_table[left_block[i]], 
                                       right_quantum_table[right_block[j]]);
        if(info != -1)
        {
            num_block++;
        }
    }

    left_index = new int[num_block];
    right_index = new int[num_block];
    physics_index = new int[num_block];
    int k = 0;
    for(int i=0;i<num_left_block;++i) for(int j=0;j<num_right_block;++j)
    {   
        info = CheckQuantumTable(site, left_quantum_table[left_block[i]], 
                                       right_quantum_table[right_block[j]]);
        if(info != -1)
        {
            // why not try 
            // left_index[k] = left_block[i];
            // right_index[k] = right_block[j];
            left_index[k] = i;
            right_index[k] = j;
            physics_index[k] = info;
            k++;
        }
    } 
}

void RealTensorSpace::DefineQuantumTable()
{
    double mean_quantum;
    int min_from_left, max_from_left, min_from_right, max_from_right, min_final, max_final;
    num_quantum_ = 1;
    num_table_ = new int[num_site_pp_];
    quantum_table_ = new int** [num_site_pp_];

    for(int i=0;i<num_site_pp_;++i)
    {
        if(i == 0)
        {
            min_final = 0;
            max_finial = 0;
            mean_quantum = 0;
        }
        else
        {
            min_from_left = 0;
            max_from_left = (num_level-1)*i;
            min_from_right = n-(num_site_-i)*(num_level-1);
            max_from_right = n;
            min_final = max(min_from_left, min_from_right);
            max_final = min(max_from_left, max_from_right);

            mean_quantum = num_boson/num_site*i;
        }
        num_tale_[i] = max_final - min_final+1;
        quantum_table_[i] = new int* [num_table_[i]];
        for(int j=0;j<num_table_[i];++j)
            quantum_table_[i][j] = new int[num_quantum_];
        for(int j=min_final;j<max_final;++j)
            quantum_table_[i][j-min_final][0] = j;
        ReorderQuantumTable(num_quantum_, num_table_[i], quantum_table_[i], &mean_quantum);
    }
}

void RealTensorSpace::ReorderQuantumTable(int num_quantum, int num_table, int** quantum_table, 
                                          double *mean_quantum)
{
    double* deviate;
    int* index;
    deviate = new double[num_table];
    
    for(int i=0;i<num_table;++i)
    {
        deviate[i] = 0.0;
        index[i] = i;
        for(int j=0;j<num_quantum;++j)
            deviate[i]+=pow(quantum_table[i][j]-mean_quantum[j],2);
    }
    QuickSort(deviate, index, 0, num_table-1, 1);
    ReorderRelevantArray(num_quantum, num_table, quantum_table, index);
}

int RealTensorSpace::CheckQuantumTable(int site, int* left_table, int* right_table)
{
    int info = -1;

    for(int i=0;i<num_level;++i)
    {
        if(left_table[0]+i == right_table[0])
        {
            info = i;
            break
        }
    }
    return info;
}



