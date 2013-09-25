#ifndef _CONVOLUTION_H_
#define _CONVOLUTION_H_

#include <math.h>
#include "complex.h"
#include "fft.h"

//
// Computes the convolution of a[0..(n-1)] and b[0..(to_fill-1)]
// and writes the result in result[0..(to_fill-1)].
// The computation is done "naively" i.e. without using FFT's.
// Note that b and result can be the same vector.
//
template< class T >
    void
    conv(
        T* a, T* b, T* result, int n, int to_fill ) {

        for( int i = to_fill - 1; i >= 0; i-- ){

            result[ i ] = a[ 0 ] * b[ i ];

            for( int j = 1; j <= MIN(i, n-1); j++ )
                result[ i ] += a[ j ] * b[ i - j ];
        }
    }

// Returns the self-convolution of a[0..(n-1)] 
// to the given power.
// The computation is done "naively".
//
template< class T >
    T*
    self_conv(
        T* a, int n, int power ) {

        T* result = new T[ n * power ];
        for( int i = 0; i < n; i++ )
            result[ i ] = a[ i ];

        for( int p = 2; p <= power; p++ )
            conv< T > (
                a, result, result, n, n * p );

        return result;
    }

//
// Computes the fft of a[0..(n-1)].
// Based on Numerical recipes code.
//
inline void
fftnr(
    complex* a, int n, int sign ) {

    fft(
        ( ( double* )a ) - 1, n, sign );
}

//
// Returns the cyclic self-convolution of a[0..(n-1)] 
// to the given power with the given period (should be a power of 2).
// Computed using FFTs and based on Numerical recipes code.
//
double*
self_conv_fftnr(
    double* a, int n, int power, int period );

//
// Returns the self-convolution of a[0..(period-1)] 
// to the all the powers between min_power and max_power, with the given period.
// Computed using FFTs and based on Numerical recipes code.
//
double**
multi_self_conv_fftnr(
    double* a, int min_power, int max_power, int period );

//
// FFTW code. To be compiled if NR = 0.
//

#if NR == 0

#include <fcntl.h>
#include <unistd.h>
#include "fftw3.h"

//
// Returns an appropriate structure for
// read and write locks to pass to the function fcntl.
//
inline struct flock*
file_lock(
    short type, short whence ) {

    static struct flock ret;
    ret.l_type = type;
    ret.l_start = 0;
    ret.l_whence = whence;
    ret.l_len = 0;
    ret.l_pid = getpid();
    return &ret;
}

//
// For importing the wisdom file for FFTW
// before computing the first fft.
//
static int wisdom_init = 0;
inline void
init_fftw() {

    if( !wisdom_init ){
        FILE* wisdom = fopen(
            "wisdom", "r" );
        int w = open(
            "wisdom", O_RDONLY );
        fcntl(
            w, F_SETLKW, file_lock(
                F_RDLCK, SEEK_SET ) );
        fftw_import_wisdom_from_file(
            wisdom );
        fcntl(
            w, F_SETLKW, file_lock(
                F_RDLCK, SEEK_SET ) );
        close(
            w );
        fclose(
            wisdom );
        wisdom_init = 1;
    }
}

//
// To update the wisdom file with newly acquired wisdom.
//
inline void
save_fftw() {

    FILE* wisdom = fopen(
        "wisdom", "w" );
    int w = open(
        "wisdom", O_WRONLY );
    fcntl(
        w, F_SETLKW, file_lock(
            F_WRLCK, SEEK_SET ) );
    fftw_export_wisdom_to_file(
        wisdom );
    fcntl(
        w, F_SETLKW, file_lock(
            F_WRLCK, SEEK_SET ) );
    close(
        w );
    fclose(
        wisdom );
}

//
// Computes the fft of a[0..(n-1)].
// Based on fftw code.
//
inline void
fftw(
    complex* a, int n, int sign ) {

    init_fftw();
    fftw_complex in[ n ];
    fftw_complex out[ n ];
    fftw_plan fft = fftw_plan_dft_1d(
        n, in, out, sign, FFTW_MEASURE );
    save_fftw();

    for( int i = 0; i < n; i++ ){
        in[ i ][ 0 ] = a[ i ].real();
        in[ i ][ 1 ] = a[ i ].imag();
    }
    fftw_execute(
        fft );

    for( int i = 0; i < n; i++ )
        a[ i ].set(
            out[ i ][ 0 ], out[ i ][ 1 ] );
}

//
// Returns the cyclic self-convolution of a[0..(n-1)] 
// to the given power with the given period.
// Computed using FFTs and based on fftw code.
//
double*
self_conv_fftw(
    double* a, int n, int power, int period );

//
// Returns the self-convolution of a[0..(period-1)] 
// to the all the powers between min_power and max_power, with the given period.
// Computed using FFTs and based on fftw code.
//
double**
multi_self_conv_fftw(
    double* a, int min_power, int max_power, int period );

#endif

#endif
