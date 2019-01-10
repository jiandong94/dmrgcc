#include "util/general.h"
#include "util/input.h"
int gettoken(istream& is, string& s);
int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        printf("Usage: %s inputfile", argv[0]);
        return 0;
    }
    auto input = InputGroup(argv[1], "input");
    auto N = input.getInt("N");
    cout << "N = " << N << endl;
    auto Jx = input.getDouble("Jx");
    cout << "Jx = " << Jx << endl;
    auto PBC = input.getBool("PBC");
    cout << "PBC = " << PBC << endl;
    auto NupPath = input.getString("NupPath");
    cout << "NupPath = " << NupPath << endl;
    
    
    N = input.getInt("N", 2);
    cout << "N = " << N << endl;
    auto N_test = input.getInt("N_test", 2);
    cout << "N_test = " << N_test << endl;
    
    auto nsweeps = input.getInt("nsweeps");
    auto table = InputGroup(input, "sweeps");
    if(!(table.GotoGroup()))
        error("Couldn't find table " + table.name());
    //string s;
    //gettoken(table.file(), s);
    //cout << s << endl;
    table.SkipLine(); // we have a table key
    int *maxm, *minm, *niter;
    double *cutoff, *noise;
    maxm = new int[nsweeps];
    minm = new int[nsweeps];
    niter = new int[nsweeps];
    cutoff = new double[nsweeps];
    noise = new double[nsweeps];
    auto nlast = nsweeps;
    for(int i=0;i<nsweeps;++i)
    {
        table.file() >> maxm[i] >> minm[i] >> cutoff[i] >> niter[i] >> noise[i];
        if(maxm[i] == 0)
        {
            nlast = i-1;
            break;
        }
    }
    for(int i=nlast+1;i<nsweeps;++i)
    {
        maxm[i] = maxm[nlast];
        minm[i] = minm[nlast];
        cutoff[i] = cutoff[nlast];
        niter[i] = niter[nlast];
        noise[i] = noise[nlast];
    }
    cout << "Sweeps:" << endl;
    for(int i=0;i<nsweeps;++i)
    {
        cout << i+1 << " " << "Maxm=" << maxm[i] << ", " <<
                              "Minm=" << minm[i] << ", " <<
                              "Cutoff=" << cutoff[i] << ", " <<
                              "Niter=" << niter[i] << ", " <<
                              "Noise=" << noise[i] << endl;
    }
    
    return 0;
}
