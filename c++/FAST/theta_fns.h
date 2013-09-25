#ifndef _THETA_OBJS_H_
#define _THETA_OBJS_H_

#include "extended_exponent.h"
#include "utils.h"

//
// Namespace for the variables and functions
// to use to compute theta in the sFFT and csFFT algorithm
//
class thetaFn {

private:

    //
    // log_pmf = Log of single column lattice pmf
    // step = Step size for the lattice
    // s = LLR score for which theta is optimized
    // Q = Lattice size
    //
    static double* log_pmf;
    static double step;
    static double s;
    static int Q;

public:

    //
    // To initialize the namespace
    //
    static void
    init(
        double* given_log_pmf, double given_step, double given_s, int given_Q ) {

        log_pmf = given_log_pmf;
        step = given_step;
        s = given_s;
        Q = given_Q;
    }

    //
    // Returns -theta*s + log(M(theta))
    //
    static double
    minimizationFn(
        double theta ) {

        double log_mgf = LOGZERO;
        for( int It = 0; It < Q; It++ )
            if( log_pmf[ It ] != LOGZERO )
                log_mgf = log_sum(
                    log_mgf, log_pmf[ It ] + theta * step * It );

        return -theta * s + log_mgf;
    }

    //
    // Returns sigma(theta)
    //
    static double
    get_sigma(
        double theta ) {

        //
        // Let p[s] = exp(log_pmf[s/step]), then,
        // log_m0 = log(\sum_s p[s]exp(theta*s))
        // log_m1 = log(\sum_s s*p[s]exp(theta*s))
        //
        double log_m0;
        double log_m1;
        log_m0 = log_m1 = LOGZERO;

        for( int i = 0; i < Q; i++ ){

            log_m0 = log_sum(
                log_m0, log_pmf[ i ] + theta * step * i );
            log_m1 = log_sum(
                log_m1, log_pmf[ i ] + theta * step * i + log(
                    step * i ) );
        }

        //
        // sigma = \sum_s (s-m1/m0)^2 p[s]exp(theta*s)/m0
        //
        double sigma = 0;
        double norm_m1 = exp(
            log_m1 - log_m0 );
        for( int i = 0; i < Q; i++ )
            sigma += pow(
                i * step - norm_m1, 2 ) * exp(
                log_pmf[ i ] + theta * i * step - log_m0 );

        return sqrt(
            sigma );
    }

    static double
    get_s() {
        return s;
    }
};

#endif
