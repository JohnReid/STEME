#ifndef _EXTENDED_EXPONENT_H_
#define _EXTENDED_EXPONENT_H_

class ee_t;

#include <math.h>
#include "utils.h"
#include <iostream>

using namespace std;

//
// This should probably change with the precision
// of the machine.
//
static const int max_shift = 60;
#define EXP1 2.71828182845905

//
// shift[i] = exp(i) for i in [0..max_shift-1]
//
inline double*
init() {

    double* shift = new double[ max_shift ];
    for( int i = 0; i < max_shift; i++ )
        shift[ i ] = exp(
            ( double )i );
    return shift;
}
static double* shift = init();

//
// Class for representing values > 0 with very small
// or very large exponents
//
class ee_t {

private:

    //
    // The value represented is base*exp(exponent)
    //
    double base;
    int exponent;

    ee_t
    normalize(
        double given_base, int given_exponent ) {

        ee_t result;
        if( given_base >= EXP1 ){
            result.base = given_base / EXP1;
            result.exponent = given_exponent + 1;
        }else{
            result.base = given_base;
            result.exponent = given_exponent;
        }
        return result;
    }

public:

    //
    // Functions for setting and getting log of the represented value.
    //
    double
    log_get() {
        return ( base == 0
            ? LOGZERO : exponent + log(
                base ) );
    }

    void
    log_set(
        double log_v ) {

        exponent = ( log_v <= LOGZERO
            ? -10000 : ( int )log_v );
        base = ( log_v <= LOGZERO
            ? 0 : exp(
                log_v - exponent ) );
    }

    ee_t() {
        base = 0;
        exponent = -10000;
    }
    ee_t(
        double log_v ) {
        log_set(
            log_v );
    }

    double
    get() {
        return base * exp(
            ( double )exponent );
    }

    ee_t&
    operator=(
        const ee_t& a ) {
        base = a.base;
        exponent = a.exponent;
        return *this;
    }

    ee_t&
    operator=(
        double log_v ) {
        log_set(
            log_v );
        return *this;
    }

    int
    operator==(
        double v ) {
        return ( v == 0
            ? ( base == 0 ) : ( get() == v ) );
    }

    int
    operator!=(
        double v ) {
        return ( v == 0
            ? ( base != 0 ) : ( get() != v ) );
    }

    ee_t
    operator*(
        const ee_t& a ) {
        return normalize(
            base * a.base, exponent + a.exponent );
    }

    ee_t
    operator+(
        const ee_t& a ) {

        //
        // Note: This special case is handled in this way
        // to optimize the runtime
        //
        if( base == 0 )
            return a;

        ee_t result;
        int diff = exponent - a.exponent;

        //
        // Case 1: a is smaller => a.base needs to be shifted
        //
        if( diff > 0 )
            if( diff < max_shift ){
                result.exponent = exponent;
                result.base = base + a.base / shift[ diff ];
                if( base >= EXP1 ){
                    result.base /= EXP1;
                    result.exponent++;
                }
            }else
                result = *this;
        //
        // Case 2: a is larger => base needs to be shifted
        //
        else if( diff <= 0 ){
            result.exponent = a.exponent;

            if( diff > -max_shift ){
                result.base = base / shift[ -diff ] + a.base;
                if( result.base >= EXP1 ){
                    result.base /= EXP1;
                    result.exponent++;
                }
            }else
                result.base = a.base;
        }

        return result;
    }

    ee_t&
    operator+=(
        const ee_t& a ) {
        *this = *this + a;
        return *this;
    }

};

inline double
log(
    ee_t& a ) {
    return a.log_get();
}

#endif
