#include "dmrg/space/bose_hubbard_space.h"

int main()
{
    cout << "=================================" << endl;
    cout << "       Test RealTensorSpace      " << endl;
    cout << "=================================" << endl;

    RealTensorSpace *tensor_space = new BoseHubbardSpace(8,1,6,2);
    tensor_space->InitializeTensorSpace(2,2);
    tensor_space->PrintTensorSpace();
    return 0;
}
