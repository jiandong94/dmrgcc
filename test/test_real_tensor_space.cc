#include "model/space/real_spinless_bose_square_space.h"

int main()
{
    cout << "=================================" << endl;
    cout << "       Test RealTensorSpace      " << endl;
    cout << "=================================" << endl;

    RealTensorSpace *tensor_space = new RealSpinlessBoseSquareSpace(8,1,6,2);
    tensor_space->DefineTensorSpace(2,2);
    tensor_space->PrintTensorSpace();
    return 0;
}
