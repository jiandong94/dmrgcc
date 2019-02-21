#ifndef TENSOR_COMPLEX_H
#define TENSOR_COMPLEX_H

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;
using std::ostream;

class Complex
{
    public:

    double real_;
    double imag_;

    Complex();

    Complex(double real, double imag);

    // double -> Complex
    Complex(double real);

    Complex(const Complex& z);

    ~Complex() {};

    friend double Real(const Complex& z);

    friend double Imag(const Complex& z);

    friend double Norm(const Complex& z);

    friend double SquareNorm(const Complex& z);

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
    friend Complex operator + (const Complex &z);
    friend Complex operator - (const Complex &z);

    friend bool operator < (const Complex& z1, const Complex& z2);
    friend bool operator <= (const Complex& z1, const Complex& z2);
    friend bool operator > (const Complex& z1, const Complex& z2);
    friend bool operator >= (const Complex& z1, const Complex& z2);
    friend bool operator == (const Complex& z1, const Complex& z2);
    friend bool operator != (const Complex& z1, const Complex& z2);
    
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

inline Complex::Complex(double real)
{
    real_ = real;
    imag_ = 0.0;
}

inline Complex::Complex(const Complex& z)
{
    real_ = z.real_;
    imag_ = z.imag_;
}

inline double Real(const Complex& z)
{
    return z.real_;
}

inline double Imag(const Complex& z)
{
    return z.imag_;
}

inline double Norm(const Complex& z)
{
    return sqrt(z.real_*z.real_ + z.imag_*z.imag_);
}

inline double SquareNorm(const Complex& z)
{
    return z.real_*z.real_ + z.imag_*z.imag_;
}

inline Complex Conj(const Complex& z)
{
    return Complex(z.real_, -z.imag_);
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

inline Complex operator+ (const Complex& z)
{
    return Complex(z.real_, z.imag_);
}

inline Complex operator- (const Complex& z)
{
    return Complex(-z.real_, -z.imag_);
}

inline bool operator < (const Complex& z1, const Complex& z2)
{
    return ((z1.real_ < z2.real_) || ((z1.real_ == z2.real_) && (z1.imag_ < z2.imag_)));
}

inline bool operator <= (const Complex& z1, const Complex& z2)
{
    return ((z1.real_ < z2.real_) || ((z1.real_ == z2.real_) && (z1.imag_ <= z2.imag_)));
}

inline bool operator > (const Complex& z1, const Complex& z2)
{
    return ((z1.real_ > z2.real_) || ((z1.real_ == z2.real_) && (z1.imag_ > z2.imag_)));
}

inline bool operator >= (const Complex& z1, const Complex& z2)
{
    return ((z1.real_ > z2.real_) || ((z1.real_ == z2.real_) && (z1.imag_ >= z2.imag_)));
}

inline bool operator == (const Complex& z1, const Complex& z2)
{
    return ((z1.real_ == z2.real_) && (z1.imag_ == z2.imag_));
}

inline bool operator != (const Complex& z1, const Complex& z2)
{
    return ((z1.real_ != z2.real_) || (z1.imag_ != z2.imag_));
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
