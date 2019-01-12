#include "dmrg/real_tensor_rundmrg.h"

RealTensorRundmrg::RealTensorRundmrg(RealTensorSpace* space, RealTensorHamiltonian* hamiltonian, 
        bool disk_cache, string cache_name, int num_sweep, InputGroup table)
{
    space_ = space;
    hamiltonian_ = hamiltonian_;
    
    num_site_ = space_->get_num_site();
    num_site_pp_ = num_site_+1;
    num_site_mm_ = num_site_-1;

    num_sweep_ = num_sweep*2;
    Sweep(table);

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
