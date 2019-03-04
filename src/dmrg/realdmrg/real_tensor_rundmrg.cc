#include "dmrg/realdmrg/real_tensor_rundmrg.h"

RealTensorRundmrg::RealTensorRundmrg(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian, 
    InputGroup& input)
{
    space_ = space;
    hamiltonian_ = hamiltonian;
  
    // cache folder
    disk_cache_ = input.getBool("disk_cache");
    string cache_name = input.getString("cache_name");
    strcpy(cache_name_, (char*) cache_name.c_str());
    MkdirCacheFolder(disk_cache_, cache_name_);
    
    // cache data
    cache_record_ = input.getBool("cache_record");
    string record_name = input.getString("record_name");
    strcpy(record_name_, (char*) record_name.c_str());
    
    // resume data
    cache_resume_ = input.getBool("cache_resume");
    string resume_name = input.getString("resume_name");
    strcpy(resume_name_, (char*) resume_name.c_str());
    
    // record process
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
    // network
    network_ = new RealTensorNetwork(space_, hamiltonian_);
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

    for(int i=0;i<num_sweep_;++i)
        noise_factor_[i] = pow(0.1, noise_factor_[i]);
}

void RealTensorRundmrg::Run()
{
    RealTensorLanczos* lanczos;
    ofstream process_file;
    double start_time, end_time;
    double energy;
    char label[512], tensor_name[512], record_name[512];
    int dummy;
    cout.setf(ios::fixed);
    cout.precision(12);
    
    start_time = get_wall_time();
    //process_file.open(process_name_, ios::binary|ios::app|ios::out);
    //process_file.precision(14);
    
    // right canonical tensor space
    // compute right contraction
    Initialize();
    lanczos = new RealTensorLanczos(space_, hamiltonian_, network_);
    
    // begin sweep
    int left = 0;
    int right = 1;
    
    for(int i=0;i<num_sweep_/2;++i)
    {
        cout<<"------------------------ sweep="<<i<<" ------------------------"<<endl;
        // sweep from left to right
        space_->DefineTensorSpaceParameter(max_block_[2*i], max_dim_[2*i], 
                canonical_precision_[2*i], noise_factor_[2*i]);
        
        for(int j=0;j<num_site_mm_;++j)
        {
            lanczos->LanczosMethod(num_iter_[2*i], j, energy);
            cout << "energy: " << energy << endl;    
            network_->ResetTensorNetwork(right, j);
            network_->RemoveTensorNetwork(right, j);

            // subspace expansion
            network_->ExpanTensorNetwork(left, j);
            space_->ExpanTensorSpace(left, j);
            
            // left canonical
            space_->CanonicalTensorSpace(left, j);
            
            // compute L tensor 
            network_->ComputeTensorNetwork(left, j);

            if(disk_cache_ == true)
            {
                space_->RecordTensorSpace(j);
                space_->ResetTensorSpace(j);

                network_->RecordTensorNetwork(left, j);
                network_->ResetTensorNetwork(left, j);

                space_->ResumeTensorSpace(j+1);
                network_->ResumeTensorNetwork(right, j+1);
            }

            space_->MatchTensorSpace(left, j+1);
            space_->MergeTensorSpace(left, j+1);
        }
         
        cout << "right energy = " << energy << endl;

        // sweep from right to left
        space_->DefineTensorSpaceParameter(max_block_[2*i+1], max_dim_[2*i+1], 
                canonical_precision_[2*i+1], noise_factor_[2*i+1]);
        
        for(int j=num_site_mm_;j>0;--j)
        {
            lanczos->LanczosMethod(num_iter_[2*i+1], j, energy);

            network_->ResetTensorNetwork(left, j);
            network_->RemoveTensorNetwork(left, j);

            // subspace expansion
            network_->ExpanTensorNetwork(right, j);
            space_->ExpanTensorSpace(right, j);
            
            // right canonical
            space_->CanonicalTensorSpace(right, j);
            
            // compute R tensor 
            network_->ComputeTensorNetwork(right, j);

            if(disk_cache_ == true)
            {
                space_->RecordTensorSpace(j);
                space_->ResetTensorSpace(j);

                network_->RecordTensorNetwork(right, j);
                network_->ResetTensorNetwork(right, j);

                space_->ResumeTensorSpace(j-1);
                network_->ResumeTensorNetwork(left, j-1);
            }

            space_->MatchTensorSpace(right, j-1);
            space_->MergeTensorSpace(right, j-1);
        }
        
        cout << "left energy = " << energy << endl;

        if(disk_cache_ == true)
        {
            sprintf(label, "rm -rf %s/TensorSpace_sweep_*.dat", cache_name_);
            dummy = system(label);
 
            sprintf(tensor_name, "%s/TensorSpace_sweep_%d.dat", cache_name_, i);
            space_->RecordTensorSpace(0);
            space_->WriteTensorSpace(tensor_name);
            space_->ResumeTensorSpace(0);
        }

        if(cache_record_ == true)
        {
            sprintf(record_name, "%s/%s", cache_name_, record_name_);
            space_->WriteTensorSpace(record_name);
        }

        if(disk_cache_ == true)
        {
            sprintf(label, "rm -rf %s/TensorSpace_sweep_*.dat", cache_name_);
            dummy = system(label);
        }

    }
    end_time = get_wall_time();
    cout << "Total Time : " << (double)(end_time-start_time) << " s" << endl;

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
