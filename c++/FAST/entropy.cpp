#include <iostream>
#include <math.h>
#include "utils.h"
#include "theta_fns.h"
#include "minimize.h"
#include "convolution.h"
#include "entropy.h"

using namespace std;

//
// Function: Computes the lattice pmf (log-values in an array of size
//           L*Q) for a sum of scores based on the sFFT algorithm.
//
// Parameters: log_pmf = Log of lattice pmf for single sample
//             L = Number of i.i.d. samples
//             s = Score for which theta is optimized
//             Q = Lattice size for single sample
//             step = Step size for the lattice
//             theta = Optimal theta shift
//             log_mgf = Log of the mgf of the pmf at theta
//
// If theta is set to NOT_COMPUTED then it is computed and otherwise
// the given theta is used.
//
double*
sFFT(
    double* log_pmf, int L, double s, int Q, double step, double& theta,
    double& log_mgf ) {

    //
    // Computing theta
    //
    if( theta == NOT_COMPUTED ){
        thetaFn::init(
            log_pmf, step, s / ( ( double )L ), Q );
        log_mgf = minimize(
            thetaFn::minimizationFn, theta );
        log_mgf += thetaFn::get_s() * theta;
        if( DEBUG )
            cout << "theta: " << theta << endl << "log(MGF): " << log_mgf
                            << endl;
    }

    //
    // Shifting the pmf
    //
    double* pmf = log_pmf;
    for( int It = 0; It < Q; It++ )
        pmf[ It ] = exp(
            log_pmf[ It ] + theta * step * It - log_mgf );

    //
    // Computing the convolution of the shifted pmf
    //
#if NR
    int period = RND_UP_TO_POWER_2(L*Q);
    double* pmfL = self_conv_fftnr(pmf, Q, L, period);
#else
    int period = L * Q;
    double* pmfL = self_conv_fftw(
        pmf, Q, L, period );
#endif

    //
    // Unshifting the result
    //
    double* log_pmfL = pmfL;
    for( int It = 0; It < period; It++ )
        log_pmfL[ It ] = ( pmfL[ It ] > 0
            ? log(
                pmfL[ It ] ) : LOGZERO ) - theta * step * It + L * log_mgf;

    return log_pmfL;
}

const int ksigma = 4;

//
// Function: Computes the lattice pmf (log-values in an array of size
//           L*Q) for a sum of scores based on the csFFT algorithm.
//
// Parameters: log_pmf = Log of lattice pmf for single sample
//             L = Number of i.i.d. samples
//             s = Score for which theta is optimized
//             Q = Lattice size for single sample
//             step = Step size for the lattice
//             theta = Optimal theta shift
//             log_mgf = Log of the mgf of the pmf at theta
//
// If theta is set to NOT_COMPUTED then it is computed and otherwise
// the given theta is used.
//
double*
csFFT(
    double* log_pmf, int L, double s, int Q, double step, double& theta,
    double& log_mgf ) {

    //
    // Computing theta
    //
    thetaFn::init(
        log_pmf, step, s / ( ( double )L ), Q );
    if( theta == NOT_COMPUTED ){
        log_mgf = minimize(
            thetaFn::minimizationFn, theta );
        log_mgf += thetaFn::get_s() * theta;
        if( DEBUG )
            cout << "theta: " << theta << endl << "log(MGF): " << log_mgf
                            << endl;
    }

    //
    // Computing the period
    //
    double sigma = thetaFn::get_sigma(
        theta );
    int Lprime = ( int )ceil(
        ksigma * sigma * sqrt(
            ( double )L ) / ( step * Q ) );
    Lprime = MIN(L, MAX(1, Lprime));

    if( DEBUG )
        cout << "sigma: " << sigma << endl << "Lprime: " << Lprime << endl;

    //
    // Boosting (Not implemented: boosting for the left tail)
    //
    if( theta > 0 && L * Q * step - s + L * 100 <= Lprime * Q * step / 2.0 ){

        theta += 20.0 / ( L * Q * step - s );
        sigma = thetaFn::get_sigma(
            theta );
        Lprime
                        = ( int )MAX(1, ceil(ksigma*sigma*sqrt((double) L)/(step*Q)));
        Lprime = MIN(L, MAX(1, Lprime));
        log_mgf = thetaFn::minimizationFn(
            theta ) + thetaFn::get_s() * theta;

        if( DEBUG )
            cout << "Boosted theta: " << theta << endl << "Boosted sigma: "
                            << sigma << endl << "Boosted log(MGF): " << log_mgf
                            << endl << "Boosted Lprime: " << Lprime << endl;
    }

    //
    // Shifting the pmf
    //
    double* pmf = log_pmf;
    for( int It = 0; It < Q; It++ )
        pmf[ It ] = exp(
            log_pmf[ It ] + theta * step * It - log_mgf );

    //
    // Computing the convolution of the shifted pmf
    //
#if NR
    int period = RND_UP_TO_POWER_2(Lprime*Q);
    double* pmfL = self_conv_fftnr(pmf, Q, L, period);
#else
    int period = Lprime * Q;
    double* pmfL = self_conv_fftw(
        pmf, Q, L, period );
#endif

    //
    // Unshifting the result
    //
    double* log_pmfL = new double[ L * Q ];
    for( int It = 0; It < L * Q; It++ )
        if( It * step > s - ( period / 2 ) * step && It * step < s + ( period
                        / 2 ) * step )
            log_pmfL[ It ] = ( pmfL[ It % period ] > 0
                ? log(
                    pmfL[ It % period ] ) : LOGZERO ) - theta * step * It + L
                            * log_mgf;
        else
            log_pmfL[ It ] = LOGZERO;

    return log_pmfL;
}

//
// Function: Computes the p-value for a sum of scores based on the lattice-refinement 
//           version of the csFFT algorithm.
//
// Parameters: log_pmf = Log of lattice pmf for single sample
//             L = Number of i.i.d. samples
//             s = Score for which theta is optimized
//             Q = Lattice size for single sample
//             step = Step size for the lattice
//             accuracy = Desired number of decimal places of accuracy for the p-value
//             theta = Optimal theta shift
//             log_mgf = Log of the mgf of the pmf at theta
//             log_max_pval = Log of an upper bound on the p-value
//             log_min_pval = Log of a lower bound on the p-value
//
// If theta is set to NOT_COMPUTED then it is computed and otherwise
// the given theta is used.
//
void
refine_csFFT(
    double* log_pmf, int L, double s, int Q, double step, double accuracy,
    double& theta, double& log_mgf, double& log_max_pval, double& log_min_pval ) {

    //
    // Computing theta
    //
    thetaFn::init(
        log_pmf, step, s / ( ( double )L ), Q );
    if( theta == NOT_COMPUTED ){
        log_mgf = minimize(
            thetaFn::minimizationFn, theta );
        log_mgf += thetaFn::get_s() * theta;
        if( DEBUG )
            cout << "theta: " << theta << endl << "log(MGF): " << log_mgf
                            << endl;
    }

    double sigma = thetaFn::get_sigma(
        theta );
    int
                    period_start =
                                    MAX(2*Q, (int) ceil(2*sigma*sqrt((double) L)/step));
    int max_multiple_of_period = RND_UP_TO_POWER_2(L*Q/period_start);

    //
    // Can switch to this version if greater efficiency is needed using FFTW
    //
    // period_start = int(L*Q/double(max_multiple_of_period));
    //
    period_start = RND_UP_TO_POWER_2(L*Q/double(max_multiple_of_period));

    if( DEBUG )
        cout << "period start: " << period_start << " " << "max multiple: "
                        << max_multiple_of_period << endl;

    //
    // Shifting the pmf
    //
    double* pmf = log_pmf;
    for( int It = 0; It < Q; It++ )
        pmf[ It ] = exp(
            log_pmf[ It ] + theta * step * It - log_mgf );

    complex* Dpmf = new complex[ period_start * max_multiple_of_period ];
    complex* Dpmf_partial = new complex[ period_start * max_multiple_of_period ];

    for( int It = 0; It < Q; It++ )
        Dpmf_partial[ It % period_start ] += complex(
            pmf[ It ] );

#if NR
    fftnr(Dpmf_partial, period_start, 1);
#else
    fftw(
        Dpmf_partial, period_start, 1 );
#endif

    for( int It = 0; It < period_start; It++ )
        Dpmf[ It ] = Dpmf_partial[ It ] ^= L;

    if( theta > 0 )
        compute_phi_pvalue(
            Dpmf, step, period_start, s, int(
                s / step + period_start / 2 ) + 1, L * 4, log_min_pval,
            log_max_pval, theta, L * log_mgf );
    else
        compute_phi_pvalue(
            Dpmf, step, period_start, s - period_start * step / 2, int(
                s / step ) + 1, L * 4, log_min_pval, log_max_pval, theta, L
                            * log_mgf );

    double old_log_pval = ( log_min_pval + log_max_pval ) / 2;
    double log_pval;

    for( int i = 2; i <= max_multiple_of_period; i *= 2 ){

        for( int It = period_start - 1; It >= 0; It-- )
            Dpmf[ 2 * It ] = Dpmf[ It ];

        double omg = pi / period_start;
        for( int It = 0; It < period_start; It++ )
            Dpmf_partial[ It ] = complex();
        for( int It = 0; It < Q; It++ )
            Dpmf_partial[ It % period_start ] += complex(
                pmf[ It ] ) * complex(
                cos(
                    It * omg ), sin(
                    It * omg ) );

#if NR
        fftnr(Dpmf_partial, period_start, 1);
#else
        fftw(
            Dpmf_partial, period_start, 1 );
#endif

        for( int It = period_start - 1; It >= 0; It-- )
            Dpmf[ 2 * It + 1 ] = Dpmf_partial[ It ] ^ L;

        period_start *= 2;

        //
        // This version needs to be debugged.
        //
        // if(theta > 0)
        //   compute_phi_pvalue(Dpmf, step, period_start, s, int(s/step+period_start/2)+1, L*4, log_min_pval, log_max_pval, theta, L*log_mgf);
        // else
        //   compute_phi_pvalue(Dpmf, step, period_start, s-period_start*step/2, int(s/step)+1, L*4, log_min_pval, log_max_pval, theta, L*log_mgf);

        complex* Dpmf_copy = new complex[ period_start ];
        for( int j = 0; j < period_start; j++ )
            Dpmf_copy[ j ] = Dpmf[ j ];

#if NR
        fftnr(Dpmf_copy, period_start, -1);
#else
        fftw(
            Dpmf_copy, period_start, -1 );
#endif

        //
        // Unshifting the result
        //
        double* log_pmfL = new double[ L * Q ];
        for( int It = 0; It < L * Q; It++ )
            if( It * step > s - ( period_start / 2 ) * step && It * step < s
                            + ( period_start / 2 ) * step )
                log_pmfL[ It ] = ( Dpmf_copy[ It % period_start ].real() > 0
                    ? log(
                        Dpmf_copy[ It % period_start ].real() ) : LOGZERO )
                                - theta * step * It + L * log_mgf - log(
                    period_start );
            else
                log_pmfL[ It ] = LOGZERO;

        compute_theta_pvalues(
            log_pmfL, step, L * Q, s, L * 4, log_min_pval, log_max_pval, theta );

        delete[] ( log_pmfL );
        delete[] ( Dpmf_copy );

        log_pval = ( log_min_pval + log_max_pval ) / 2;

#if DEBUG
        cout << i << ": " << old_log_pval << " " << log_pval << " " << log_min_pval << " " << log_max_pval << endl;
#endif

        if( -log(
            fabs(
                log_pval - old_log_pval ) ) / log(
            10 ) > accuracy )
            break;

        old_log_pval = log_pval;
    }

#if DEBUG
    cout << "period end: " << period_start << endl;
#endif

    //
    // For the "compute_phi_pvalue" based version"
    //
    // if(theta < 0) {
    //
    // double log_min_pval_copy = log_min_pval;
    // log_min_pval = log(1-exp(log_max_pval));
    // log_max_pval = log(1-exp(log_min_pval_copy));
    // }
    //

    return;
}
