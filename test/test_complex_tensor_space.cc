#include "model/space/complex_spinless_bose_square_space.h"

int main()
{
    cout << "=================================" << endl;
    cout << "    Test ComplexTensorSpace      " << endl;
    cout << "=================================" << endl;

    ComplexTensorSpace *tensor_space = new ComplexSpinlessBoseSquareSpace(8,1,6,2);
    tensor_space->DefineTensorSpace(2,2);
    tensor_space->PrintTensorSpace();
    return 0;
}
