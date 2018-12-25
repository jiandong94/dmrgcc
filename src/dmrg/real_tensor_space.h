#ifndef DMRGCC_DMRG_REAL_TENSOR_SPACE_H_
#define DMRGCC_DMRG_REAL_TENSOR_SPACE_H_

#include "dmrg/real_tensor_lattice.h"

class RealTensorSpace
{
    protected:

    // disk cache
    bool disk_cache_;
    char cache_name_[512];

    // number of sites
    int num_site_;
    int num_site_pp_;
    int num_site_mm_;

    // quantum table
    int num_quantum_;
    int *num_table_;
    int ***quantum_table_;

    // tensor lattice
    RealTensorLattice **tensor_lattice_;

    // maximun number of block
    int max_block_;
    // maximun bond dimension
    int max_dim_;
    // canonical precision
    double canonical_precision_;
    // noise factor
    double noise_factor_;

    // expan lattice
    RealTensorLattice* expan_tensor_lattice_;

    public:
    
    // destructor
    //
    virtual ~RealTensorSpace();

    // 
    //
    RealTensorLattice* get_tensor_lattice(int site);

    //
    //
    void DefineTensorSpaceCache(bool disk_cache, char *cache_name);

    //
    //
    void DefineTensorSpaceParameter(int max_block, int max_dim, 
            double canonical_precision, double noise_factor);

    
    void InitializeTensorSpace(int initial_block, int initial_dim);
    //
    //
    void PrintTensorSpace();

    //
    //
    void WriteTensorSpace(const char* tensor_space_name);

    //
    //
    void WriteTensorSpace(ofstream &tensor_space_file);

    //
    //
    void ReadTensorSpace(const char* tensor_space_name);

    //
    //
    void ReadTensorSpace(ifstream &tensor_space_file);

    //
    //
    void RecordTensorSpace(int site);

    //
    //
    void ResumeTensorSpace(int site);

    //
    //
    void ResetTensorSpace(int site);

};


#endif // DMRGCC_DMRG_REAL_TENSOR_SPACE_H_
