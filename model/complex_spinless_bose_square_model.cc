#include "dmrgcc.h"

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        printf("Usage: %s inputfile", argv[0]);
        return 0;
    }

    auto input = InputGroup(argv[1], "input");
    int num_site_x = input.getInt("num_site_x");
    int num_site_y = input.getInt("num_site_y");
    int num_boson = input.getInt("num_boson");
    int physics_dim = input.getInt("physics_dim", 2);
    double flux = input.getDouble("flux");
    double flux_value = flux*PI;
    double Jx = input.getDouble("Jx");
    double Jy = input.getDouble("Jy");
    double U = input.getDouble("U");
    double mu = input.getDouble("mu");
    bool period_x = input.getBool("period_x");
    bool period_y = input.getBool("period_y");

    // input file should also include
    // disk_cache, cache_name, cache_resume, resume_name
    // record_process, process_name, first_block, first_dim
    // num_sweep, sweep table

    ComplexTensorSpace* space;
    ComplexTensorHamiltonian* hamiltonian;
    ComplexTensorRundmrg* rundmrg;

    space = new ComplexSpinlessBoseSquareSpace(num_site_x, num_site_y, num_boson, physics_dim);
    hamiltonian = new ComplexSpinlessBoseSquareHamiltonian(num_site_x, num_site_y, physics_dim, flux_value, Jx, Jy, 
                                             U, mu, period_x, period_y);
    rundmrg = new ComplexTensorRundmrg(space, hamiltonian, input);
    rundmrg->Run();
    
    delete space;
    delete hamiltonian;
    delete rundmrg;
}
