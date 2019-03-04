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
    int physics_dim = input.getInt("physics_dim", 2);
    bool period_x = input.getBool("period_x");
    bool period_y = input.getBool("period_y");

    int num_boson[2];
    double flux[2];
    double Jx[2];
    double Jy[2];
    double U[3];
    double mu[2];

    num_boson[0] = input.getInt("Nup");
    num_boson[1] = input.getInt("Ndown");
    for(int t=0;t<2;++t)
    {
        flux[t] = input.getDouble("flux");
        Jx[t] = input.getDouble("Jx");
        Jy[t] = input.getDouble("Jy");
        mu[t] = input.getDouble("mu");
    }
    U[0] = input.getDouble("U_uu");
    U[1] = input.getDouble("U_dd");
    U[2] = input.getDouble("U_ud");

    // input file should also include
    // disk_cache, cache_name, cache_resume, resume_name
    // record_process, process_name, first_block, first_dim
    // num_sweep, sweep table

    ComplexTensorSpace* space;
    ComplexTensorHamiltonian* hamiltonian;
    ComplexTensorRundmrg* rundmrg;

    space = new ComplexSpinfulBoseSquareSpace(num_site_x, num_site_y, num_boson, physics_dim);
    hamiltonian = new ComplexSpinfulBoseSquareHamiltonian(num_site_x, num_site_y, physics_dim, flux, Jx, Jy, 
                                             U, mu, period_x, period_y);
    rundmrg = new ComplexTensorRundmrg(space, hamiltonian, input);
    rundmrg->Run();
    
    delete space;
    delete hamiltonian;
    delete rundmrg;
}
