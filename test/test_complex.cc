#include "tensor/complex.h"

int main()
{
    
    Complex complex1, complex2;
    complex1 = Complex(1.0, 1.0);
    complex2 = Complex(1.0, 2.0);
    cout << "complex1: " << complex1 << endl;
    cout << "complex2: " << complex2 << endl;

    cout << "------------------------" << endl;
    cout << "complex1.real: " << Real(complex1) << endl;
    cout << "complex1.imag: " << Imag(complex1) << endl;
    cout << "complex1.norm: " << Norm(complex1) << endl;
    cout << "complex1.squre_norm: " << SquareNorm(complex1) << endl;
    cout << "complex1.conj: " << Conj(complex1) << endl;

    cout << "-----------------------" << endl;
    cout << "complex1+complex2: " << complex1+complex2 << endl;
    cout << "complex1-complex2: " << complex1-complex2 << endl;
    cout << "complex1*complex2: " << complex1*complex2 << endl;
    cout << "complex1/complex2: " << complex1/complex2 << endl;

    cout << "-----------------------" << endl;
    complex1 = Complex(1.0, 1.0);
    complex1 += complex2;  
    cout << "complex1+=complex2 " << complex1 << endl;
    complex1 = Complex(1.0, 1.0);
    complex1 -= complex2;  
    cout << "complex1-=complex2 " << complex1 << endl;
    complex1 = Complex(1.0, 1.0);
    complex1 *= complex2;  
    cout << "complex1*=complex2 " << complex1 << endl;
    complex1 = Complex(1.0, 1.0);
    complex1 /= complex2;  
    cout << "complex1/=complex2 " << complex1 << endl;

    cout << "-----------------------" << endl;
    complex1 = Complex(1);
    complex1 += 1.0;
    cout << "complex1+=1.0 " << complex1 << endl;
    complex1 -= 1.0;
    cout << "complex1-=1.0 " << complex1 << endl;
    complex1 *= 1.0;
    cout << "complex1*=1.0 " << complex1 << endl;
    complex1 /= 1.0;
    cout << "complex1/=1.0 " << complex1 << endl;

    cout << "-----------------------" << endl;
    cout << "complex1 + 1.0" << complex1+1.0 << endl;
    cout << "complex1 - 1.0" << complex1-1.0 << endl;
    cout << "complex1 * 2.0" << complex1*2.0 << endl;
    cout << "complex1 / 2.0" << complex1/2 << endl;
    cout << "1.0 + complex1" << 1.0+complex1 << endl;
    cout << "1.0 - complex1" << 1.0-complex1 << endl;
    cout << "2.0 * complex1" << 2.0*complex1 << endl;
    cout << "2.0 / complex1" << 2.0/complex1 << endl;
    cout << "complex1 + (1.0,1.0)" << complex1+Complex(1.0,1.0) << endl;
    cout << "complex1 - (1.0,1.0)" << complex1-Complex(1.0,1.0) << endl;
    cout << "complex1 * (1.0,0.0)" << complex1*Complex(1.0,0.0) << endl;
    cout << "complex1 / (1.0,0.0)" << complex1/Complex(1.0,0.0) << endl;
    cout << "+complex" << +complex1 << endl;
    cout << "-complex" << -complex1 << endl;




    
    return 0;
}
