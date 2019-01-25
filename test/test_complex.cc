#include "tensor/complex.h"

int main()
{
    
    Complex complex1, complex2;
    complex1 = Complex(1, 1);
    complex2 = Complex(1, 2);
    cout << "complex1: " << complex1 << endl;
    cout << "complex2: " << complex2 << endl;

    cout << "------------------------" << endl;
    cout << "complex1.real: " << complex1.Real() << endl;
    cout << "complex1.imag: " << complex1.Imag() << endl;
    cout << "complex1.norm: " << complex1.Norm() << endl;
    cout << "complex1.squre_norm: " << complex1.SquareNorm() << endl;

    cout << "-----------------------" << endl;
    cout << "complex1+complex2: " << complex1+complex2 << endl;
    cout << "complex1-complex2: " << complex1-complex2 << endl;
    cout << "complex1*complex2: " << complex1*complex2 << endl;
    cout << "complex1/complex2: " << complex1/complex2 << endl;

    cout << "-----------------------" << endl;
    complex1 = Complex(1, 1);
    complex1 += complex2;  
    cout << "complex1+=complex2 " << complex1 << endl;
    complex1 = Complex(1, 1);
    complex1 -= complex2;  
    cout << "complex1-=complex2 " << complex1 << endl;
    complex1 = Complex(1, 1);
    complex1 *= complex2;  
    cout << "complex1*=complex2 " << complex1 << endl;
    complex1 = Complex(1, 1);
    complex1 /= complex2;  
    cout << "complex1/=complex2 " << complex1 << endl;
    
    return 0;
}
