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

    bool disk_cache = input.getBool("disk_cache");
    bool disk_record = input.getBool("disk_record");
    string cache_name = input.getString("cache_name");
    bool disk_resume = input.getBool("disk_resume");
    string resume_name = input.getString("resume_name");

    int first_block = input.getInt("first_block");
    int first_dim = input.getInt("first_dim");
    int num_sweep =  input.getInt("num_sweep");

    auto table = InputGroup(input, "sweeps");
    
    RealTensorSpace* space;
    RealTensorHamiltonian* hamiltonian;
    RealTensorRundmrg* rundmrg;

    space = new BoseHubbardSpace(num_site_x, num_site_y, num_boson, physics_dim);
    hamiltonian = new BoseHubbardHamiltonian(num_site_x, num_site_y, flux_value, Jx, Jy, 
                                             U, mu, period_x, period_y);
    rundmrg = new RealTensorRundmrg(space, hamiltonian, disk_cache, cache_name,
            num_sweep, table);
    
    delete space;
    delete hamiltonian;
    delete rundmrg;
}
