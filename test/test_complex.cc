#include "tensor/complex.h"

int main()
{
    
    Complex complex1, complex2;
    complex1 = Complex(1, 1);
    complex2 = Complex(1, 2);
    cout << "complex1: " << complex1 << endl;
    cout << "complex2: " << complex2 << endl;

    cout << "------------------------" << endl;
    cout << "complex1.real: " << Real(complex1) << endl;
    cout << "complex1.imag: " << Imag(complex1) << endl;
    cout << "complex1.norm: " << Norm(complex1) << endl;
    cout << "complex1.squre_norm: " << SquareNorm(complex1) << endl;

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
