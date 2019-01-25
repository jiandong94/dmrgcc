#ifndef TENSOR_COMPLEX_H
#define TENSOR_COMPLEX_H

#include "util/general.h"

class Complex
{
    public:

    double real_;
    double imag_;

    Complex();

    Complex(double real, double imag);

    Complex(const Complex& z);

    ~Complex() {};

    double Real();

    double Imag();

    double Norm();

    double SquareNorm();

    friend Complex Conj(const Complex& z);

    Complex& operator += (const Complex& z);
    Complex& operator -= (const Complex& z);
    Complex& operator *= (const Complex& z);
    Complex& operator /= (const Complex& z);
    Complex& operator += (double r);
    Complex& operator -= (double r);
    Complex& operator *= (double r);
    Complex& operator /= (double r);
    
    friend Complex operator + (const Complex &z1, const Complex& z2);
    friend Complex operator - (const Complex &z1, const Complex& z2);
    friend Complex operator * (const Complex &z1, const Complex& z2);
    friend Complex operator / (const Complex &z1, const Complex& z2);
    friend Complex operator + (const Complex &z, double r);
    friend Complex operator - (const Complex &z, double r);
    friend Complex operator * (const Complex &z, double r);
    friend Complex operator / (const Complex &z, double r);
    friend Complex operator + (double r, const Complex &z);
    friend Complex operator - (double r, const Complex &z);
    friend Complex operator * (double r, const Complex &z);
    friend Complex operator / (double r, const Complex &z);

    Complex& operator = (const Complex& z);
    Complex& operator = (double r);

    friend ostream& operator<< (ostream& str, const Complex& z);
};

inline Complex::Complex()
{
    real_ = 0.0;
    imag_ = 0.0;
}

inline Complex::Complex(double real, double imag)
{
    real_ = real;
    imag_ = imag;
}

inline Complex::Complex(const Complex& z)
{
    real_ = z.real_;
    imag_ = z.imag_;
}

inline double Complex::Real()
{
    return real_;
}

inline double Complex::Imag()
{
    return imag_;
}

inline double Complex::Norm()
{
    return sqrt(real_*real_ + imag_*imag_);
}

inline double Complex::SquareNorm()
{
    return real_*real_ + imag_*imag_;
}

inline Complex(const Complex& z)
{
    return Complex(z.real_, -z.imag);
}

//----------------------operator-------------------------
//
inline Complex& Complex::operator+= (const Complex& z)
{
    real_ += z.real_;
    imag_ += z.imag_;
    return *(this);
}

inline Complex& Complex::operator-= (const Complex& z)
{
    real_ -= z.real_;
    imag_ -= z.imag_;
    return *(this);
}

inline Complex& Complex::operator*= (const Complex& z)
{
    double real = real_*z.real_ - imag_*z.imag_;
    imag_ = real_*z.imag_ + imag_*z.real_;
    real_ = real;
    return *(this);
}

inline Complex& Complex::operator/= (const Complex& z)
{
    double k = 1.0/(z.real_*z.real_+z.imag_*z.imag_);
    double real = (real_*z.real_ + imag_*z.imag_)*k;
    imag_ = (imag_*z.real_ - real_*z.imag_)*k;
    real_ = real;
    return *(this);
}

inline Complex& Complex::operator+= (double r)
{
    real_ += r;
    return *(this);
}

inline Complex& Complex::operator-= (double r)
{
    real_ -= r;
    return *(this);
}

inline Complex& Complex::operator*= (double r)
{
    real_ *= r;
    imag_ *= r;
    return *(this);
}

inline Complex& Complex::operator/= (double r)
{
    real_ /= r;
    imag_ /= r;
    return *(this);
}

inline Complex operator+ (const Complex& z1, const Complex& z2)
{
    return Complex(z1.real_+z2.real_, z1.imag_+z2.imag_);
}

inline Complex operator- (const Complex& z1, const Complex& z2)
{
    return Complex(z1.real_-z2.real_, z1.imag_-z2.imag_);
}

inline Complex operator* (const Complex& z1, const Complex& z2)
{
    return Complex(z1.real_*z2.real_ - z1.imag_*z2.imag_, 
                   z1.real_*z2.imag_ + z1.imag_*z2.real_);
}

inline Complex operator/ (const Complex& z1, const Complex& z2)
{
    double k = 1.0/(z2.real_*z2.real_ + z2.imag_*z2.imag_);
    return Complex((z1.real_*z2.real_ + z1.imag_*z2.imag_)*k, 
                   (z1.imag_*z2.real_ - z1.real_*z2.imag_)*k);
}

inline Complex operator+ (const Complex& z, double r)
{
    return Complex(z.real_+r, z.imag_);
}

inline Complex operator- (const Complex& z, double r)
{
    return Complex(z.real_-r, z.imag_);
}

inline Complex operator* (const Complex& z, double r)
{
    return Complex(z.real_*r, z.imag_*r);
}

inline Complex operator/ (const Complex& z, double r)
{
    return Complex(z.real_/r, z.imag_/r);
}

inline Complex operator+ (double r, const Complex& z)
{
    return Complex(z.real_+r, z.imag_);
}

inline Complex operator- (double r, const Complex& z)
{
    return Complex(r-z.real_, z.imag_);
}

inline Complex operator* (double r, const Complex& z)
{
    return Complex(z.real_*r, z.imag_*r);
}

inline Complex operator/ (double r, const Complex& z)
{
    double k = 1.0/(z.real_*z.real_ + z.imag_*z.imag_);
    return Complex((r*z.real_)*k, (-r*z.imag_)*k);
}

inline Complex& Complex::operator= (const Complex& z)
{
    real_ = z.real_;
    imag_ = z.imag_;
    return *(this);
}

inline Complex& Complex::operator= (double real)
{
    real_ = real;
    imag_ = 0.0;
    return *(this);
}


inline ostream& operator<< (ostream& str, const Complex& z)
{
    str << "(" << z.real_ << "," << z.imag_ << ")";
    return str;
}
#endif //TENSOR_COMPLEX_H
