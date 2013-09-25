#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include "complex.h"
#include "fft.h"
#include "fftw3.h"
#include "convolution.h"

#include <stdio.h>

//
// Sets r to c^power 
// (for r and c, 0th entry = real part, 1st entry = imag part)
// 
void
complex_pow(
    double *c, int power, double *r ) {

    double real = 1, imag = 0, re;
    int i;

    for( i = 0; i < power; i++ ){
        re = real * c[ 0 ] - imag * c[ 1 ];
        imag = imag * c[ 0 ] + real * c[ 1 ];
        real = re;
    }
    r[ 0 ] = real;
    r[ 1 ] = imag;
    return;
}

//
// Returns the cyclic self-convolution of a[0..(n-1)] 
// to the given power with the given period (should be a power of 2).
// Computed using FFTs and based on Numerical recipes code.
//
double*
self_conv_fftnr(
    double* a, int n, int power, int period ) {

    double* result = new double[ period ];
    for( int i = 0; i < n; i++ )
        result[ i % period ] += a[ i ];

    //
    // The argument is result-1 because of
    // the peculiarities of NR code.
    //
    fft_real(
        result - 1, period, 1 );

    //
    // The result from NR code for real FFT
    // has a particular way of storing the result
    // and hence this code for computing the power of
    // the entries.
    //
    for( int i = 0; i <= 1; i++ )
        result[ i ] = pow(
            result[ i ], power );
    for( int i = 2; i < period; i += 2 )
        complex_pow(
            result + i, power, result + i );

    fft_real(
        result - 1, period, -1 );
    for( int i = 0; i < period; i++ )
        result[ i ] *= ( double )2 / period;

    return result;
}

//
// Returns the self-convolution of a[0..(period-1)] 
// to the all the powers between min_power and max_power, with the given period.
// Computed using FFTs and based on Numerical recipes code.
//
double**
multi_self_conv_fftnr(
    double* a, int min_power, int max_power, int period ) {

    //
    // fft_a = fft(a) and fft_a_power = fft_a to the current power
    //
    double* fft_a = new double[ period ];
    double* fft_a_power = new double[ period ];

    for( int i = 0; i < period; i++ )
        fft_a[ i ] = a[ i ];

    fft_real(
        fft_a - 1, period, 1 );

    for( int i = 0; i <= 1; i++ )
        fft_a_power[ i ] = pow(
            fft_a[ i ], min_power );
    for( int i = 2; i < period; i += 2 )
        complex_pow(
            fft_a + i, min_power, fft_a_power + i );

    double** result = new double*[ max_power - min_power + 1 ];

    for( int power = min_power; power <= max_power; power++ ){

        //
        // computing result[power-min_power]
        //
        result[ power - min_power ] = new double[ period ];
        for( int i = 0; i < period; i++ )
            result[ power - min_power ][ i ] = fft_a_power[ i ];

        fft_real(
            result[ power - min_power ] - 1, period, -1 );

        for( int i = 0; i < period; i++ )
            result[ power - min_power ][ i ] *= 2 / double(
                period );

        //
        // Incrementing the power for fft_a_power.
        //
        for( int i = 0; i <= 1; i++ )
            fft_a_power[ i ] *= fft_a[ i ];
        for( int i = 2; i < period; i += 2 ){
            double temp = fft_a_power[ i ] * fft_a[ i ] - fft_a_power[ i + 1 ]
                            * fft_a[ i + 1 ];
            fft_a_power[ i + 1 ] = fft_a_power[ i + 1 ] * fft_a[ i ]
                            + fft_a_power[ i ] * fft_a[ i + 1 ];
            fft_a_power[ i ] = temp;
        }
    }

    return result;
}

//
// FFTW code. To be compiled if NR = 0.
//

#if NR == 0

//
// Returns the cyclic self-convolution of a[0..(n-1)] 
// to the given power with the given period.
// Computed using FFTs and based on fftw code.
//
double*
self_conv_fftw(
    double* a, int n, int power, int period ) {

    //
    // Loading FFTW wisdom, generating plans and saving
    // new wisdom.
    //
    init_fftw();
    double* in = new double[ period ];
    fftw_complex* out = new fftw_complex[ period ];
    fftw_plan fft = fftw_plan_dft_r2c_1d(
        period, in, out, FFTW_MEASURE );
    fftw_plan ifft = fftw_plan_dft_c2r_1d(
        period, out, in, FFTW_MEASURE );
    save_fftw();

    //
    // Doing the cyclic self-convolution
    //

    for( int i = 0; i < n; i++ )
        in[ i % period ] += a[ i ];

    fftw_execute(
        fft );

    for( int i = 0; i < period; i++ )
        complex_pow(
            out[ i ], power, out[ i ] );

    fftw_execute(
        ifft );
    delete[] ( out );

    for( int i = 0; i < period; i++ )
        in[ i ] /= ( double )period;

    return in;
}

//
// Returns the self-convolution of a[0..(period-1)] 
// to the all the powers between min_power and max_power, with the given period.
// Computed using FFTs and based on fftw code.
//
double**
multi_self_conv_fftw(
    double* a, int min_power, int max_power, int period ) {

    static int curr_period = 0;
    static double* in1 = 0;
    static fftw_complex* out1 = 0;
    static fftw_complex* in2 = 0;
    static double* out2 = 0;
    static fftw_complex* out_exp = 0;

    static fftw_plan fft = 0;
    static fftw_plan ifft = 0;

    //
    // Loading FFTW wisdom, generating plans and saving
    // new wisdom.
    //

    if( period != curr_period ){

        init_fftw();
        if( curr_period != 0 ){

            delete[] ( in1 );
            delete[] ( out1 );
            delete[] ( in2 );
            delete[] ( out2 );
            delete[] ( out_exp );
        }

        in1 = new double[ period ];
        out1 = new fftw_complex[ period ];
        in2 = new fftw_complex[ period ];
        out2 = new double[ period ];
        out_exp = new fftw_complex[ period ];

        fft = fftw_plan_dft_r2c_1d(
            period, in1, out1, FFTW_MEASURE );
        ifft = fftw_plan_dft_c2r_1d(
            period, in2, out2, FFTW_MEASURE );
        save_fftw();
        curr_period = period;
    }

    //
    // out1 = fft(a) and out_exp = out1 raised to the appropriate power
    //

    for( int i = 0; i < period; i++ )
        in1[ i ] = a[ i ];

    fftw_execute(
        fft );

    for( int i = 0; i < period; i++ )
        complex_pow(
            out1[ i ], min_power, out_exp[ i ] );

    double** result = new double*[ max_power - min_power + 1 ];
    for( int power = min_power; power <= max_power; power++ ){

        //
        // computing result[power-min_power]
        //
        for( int i = 0; i < period; i++ ){
            in2[ i ][ 0 ] = out_exp[ i ][ 0 ];
            in2[ i ][ 1 ] = out_exp[ i ][ 1 ];
        }

        fftw_execute(
            ifft );

        result[ power - min_power ] = new double[ period ];
        for( int i = 0; i < period; i++ )
            result[ power - min_power ][ i ] = out2[ i ] / period;

        //
        // Incrementing the power for out_exp.
        //
        for( int i = 0; i < period; i++ ){
            double temp = out_exp[ i ][ 0 ] * out1[ i ][ 0 ] - out_exp[ i ][ 1 ]
                            * out1[ i ][ 1 ];
            out_exp[ i ][ 1 ] = out_exp[ i ][ 1 ] * out1[ i ][ 0 ]
                            + out_exp[ i ][ 0 ] * out1[ i ][ 1 ];
            out_exp[ i ][ 0 ] = temp;
        }
    }

    return result;
}

#endif
