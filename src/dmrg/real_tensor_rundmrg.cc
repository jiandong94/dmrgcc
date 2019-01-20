#include "dmrg/real_tensor_rundmrg.h"

RealTensorRundmrg::RealTensorRundmrg(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian, 
    InputGroup& input)
{
    space_ = space;
    hamiltonian_ = hamiltonian;
   
    disk_cache_ = input.getBool("disk_cache");
    string cache_name = input.getString("cache_name");
    strcpy(cache_name_, (char*) cache_name.c_str());
    MkdirCacheFolder(disk_cache_, cache_name_);
    
    cache_resume_ = input.getBool("cache_resume");
    string resume_name = input.getString("resume_name");
    strcpy(resume_name_, (char*) resume_name.c_str());
    
    record_process_ = input.getBool("record_process");
    string process_name = input.getString("process_name");
    strcpy(process_name_, (char*) process_name.c_str());

    num_site_ = space_->get_num_site();
    num_site_pp_ = num_site_+1;
    num_site_mm_ = num_site_-1;


    int first_block = input.getInt("first_block");
    int first_dim = input.getInt("first_dim");
    int num_sweep =  input.getInt("num_sweep");
    num_sweep_ = num_sweep*2;
    auto table = InputGroup(input, "sweeps");
    Sweep(table);

    // space
    space->DefineTensorSpaceCache(disk_cache_, cache_name_);
    if(cache_resume_ == true)
    {
        space_->ReadTensorSpace(resume_name_);
    }
    else
    {   
        // left canonical
        space_->DefineTensorSpace(first_block, first_dim);
    }
    
    // hamiltonian
    hamiltonian_->DefineTensorHamiltonian();
    //hamiltonian_->PrintTensorHamiltonian();
    // network
    network_ = new RealTensorNetwork(space_, hamiltonian_);
    //network_->PrintTensorNetwork();
}

RealTensorRundmrg::~RealTensorRundmrg()
{
    delete network_;

    delete[] max_dim_;
    delete[] max_block_;
    delete[] num_iter_;
    delete[] canonical_precision_;
    delete[] noise_factor_;
}

void RealTensorRundmrg::Sweep(InputGroup &table)
{
    if(!(table.GotoGroup()))
        error("Couldn't find table " + table.name());
    
    table.SkipLine(); // we have a table key
    max_dim_ = new int[num_sweep_];
    max_block_ = new int[num_sweep_];
    num_iter_ = new int[num_sweep_];
    canonical_precision_ = new double[num_sweep_];
    noise_factor_ = new double[num_sweep_];
    auto nlast = num_sweep_;
    for(int i=0;i<num_sweep_;++i)
    {
        table.file() >> max_dim_[i] >> max_block_[i] >> num_iter_[i] >> 
                        canonical_precision_[i] >> noise_factor_[i];
        if(max_dim_[i] == 0)
        {
            nlast = i-1;
            break;
        }
    }
    for(int i=nlast+1;i<num_sweep_;++i)
    {
        max_dim_[i] = max_block_[nlast];
        max_block_[i] = max_block_[nlast];
        num_iter_[i] = num_iter_[nlast];
        canonical_precision_[i] = canonical_precision_[nlast];
        noise_factor_[i] = noise_factor_[nlast];
    }
    cout << "Sweeps:" << endl;
    for(int i=0;i<num_sweep_;++i)
    {
        cout << i+1 << " " << "Maxm=" << max_dim_[i] << ", " <<
                                "Maxb=" << max_block_[i] << ", " <<
                                "Niter=" << num_iter_[i] << ", " <<
                                "Canon=" << canonical_precision_[i] << ", " <<
                                "Noise=" << noise_factor_[i] << endl;
    }
}

void RealTensorRundmrg::Run()
{
    ofstream process_file;
    double start_time, end_time;

    //process_file.open(process_name_, ios::binary|ios::app|ios::out);
    //process_file.precision(14);
    
    // right canonical tensor space
    // compute right contraction
    Initialize();
    

}

void RealTensorRundmrg::Initialize()
{
    int left, right;
    left = 0;
    right = 1;

    // initialize L and R tensors
    // record the first and the last L and R tensors
    network_->DefineTensorNetwork();
    
    // resume the last lattice
    // resume the fiest and the last L and R tensor
    if(disk_cache_ == true)
    {
        space_->ResumeTensorSpace(num_site_mm_);

        network_->ResumeTensorNetwork(right, num_site_mm_);
        network_->ResumeTensorNetwork(left, 0);
    }
    
    // right canonical lattice and compute R tensors
    // and record them
    for(int i=num_site_mm_;i>0;--i)
    {
        space_->CanonicalTensorSpace(right, i);
        network_->ComputeTensorNetwork(right, i);

        if(disk_cache_ == true)
        {
            space_->RecordTensorSpace(i);
            space_->ResetTensorSpace(i);

            network_->RecordTensorNetwork(right, i);
            network_->ResetTensorNetwork(right, i);

            space_->ResumeTensorSpace(i-1);
        }

        space_->MergeTensorSpace(right, i-1);
    }

    space_->CanonicalTensorSpace(right, 0);
}
