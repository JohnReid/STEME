#ifndef _COMPLEX_H_
#define _COMPLEX_H_

class complex;

#include <math.h>
#include "utils.h"

//
// Class for representing complex values
//
class complex {

private:

    double re;
    double im;

public:

    inline void
    set(
        double given_re, double given_im ) {
        re = given_re;
        im = given_im;
    }

    //
    // Constructors
    //
    inline
    complex() {
        set(
            0, 0 );
    }
    inline
    complex(
        double given_re ) {
        set(
            given_re, 0 );
    }
    inline
    complex(
        double given_re, double given_im ) {
        set(
            given_re, given_im );
    }
    inline
    complex(
        const complex& c ) {
        set(
            c.re, c.im );
    }

    inline complex&
    operator=(
        const complex& c ) {
        set(
            c.re, c.im );
        return *this;
    }

    inline double
    real() {
        return re;
    }
    inline double
    imag() {
        return im;
    }

    //
    // Returns the conjugate of the number
    //
    inline complex
    conj() const {
        return complex(
            re, -im );
    }

    inline complex
    operator-(
        const complex& c ) const {
        return complex(
            re - c.re, im - c.im );
    }
    inline complex
    operator-(
        double d ) const {
        return complex(
            re - d, im );
    }
    inline complex
    operator+(
        const complex& c ) const {
        return complex(
            re + c.re, im + c.im );
    }
    inline complex&
    operator+=(
        const complex& c ) {
        set(
            re + c.re, im + c.im );
        return *this;
    }

    inline complex
    operator*(
        double given_re ) const {
        return complex(
            re * given_re, im * given_re );
    }
    inline complex
    operator*(
        const complex& c ) const {
        return complex(
            re * c.re - im * c.im, re * c.im + im * c.re );
    }
    inline complex&
    operator*=(
        const complex& c ) {
        set(
            re * c.re - im * c.im, re * c.im + im * c.re );
        return *this;
    }

    inline complex
    operator/(
        const complex& c ) const {
        return complex(
            ( re * c.re + im * c.im ) / ( c.re * c.re + c.im * c.im ), ( im
                            * c.re - re * c.im ) / ( c.re * c.re + c.im * c.im ) );
    }

    inline complex
    operator/(
        double given_re ) const {
        return complex(
            re / given_re, im / given_re );
    }

    //
    // Returns (re + im*i)*exp(log_shift)
    // (computes it carefully to preserve accuracy)
    //
    inline complex
    shift(
        double log_shift ) {
        return complex(
            sign(
                re ) * exp(
                log(
                    fabs(
                        re ) ) + log_shift ), sign(
                im ) * exp(
                log(
                    fabs(
                        im ) ) + log_shift ) );
    }

    inline complex
    operator^(
        int power ) {

        complex result(
            1 );
        for( int i = 1; i <= power; i++ )
            result *= *this;

        return result;
    }

    inline complex&
    operator^=(
        int power ) {
        *this = *this ^ power;
        return *this;
    }

};

inline double
abs(
    complex a ) {
    return sqrt(
        a.real() * a.real() + a.imag() * a.imag() );
}
inline complex
exp(
    complex a ) {
    return complex(
        cos(
            a.imag() ), sin(
            a.imag() ) ) * exp(
        a.real() );
}
inline complex
sinh(
    complex a ) {
    return complex(
        sinh(
            a.real() ) * cos(
            a.imag() ), cosh(
            a.real() ) * sin(
            a.imag() ) );
}

#endif
