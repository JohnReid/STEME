#include "utils.h"
#include "HS.h"
#include "entropy.h"
#include "entropy-combined.h"

//
// Computing the p-value for the entropy score using the sFFT algorithm.
// 
entropy_return_value
HS_sFFT(
    int N, int K, int L, double s, double* pu, int Q ) {

    //
    // rv = Value to be returned.
    //
    entropy_return_value rv;
    rv.Q = Q;

    //
    // Computing the pmf for a single column.
    //
    int* index;
    int size;
    double* log_pmf = HSlist(
        N, K, pu, Q, rv.step, index, size );

    //
    // Shifted L-fold convolution.
    //
    rv.theta = NOT_COMPUTED;
    rv.log_pmfL = sFFT(
        log_pmf, L, s, Q, rv.step, rv.theta, rv.log_mgf );
    compute_theta_pvalues(
        rv.log_pmfL, rv.step, L * Q, s, L * K, rv.log_min_pvalue,
        rv.log_max_pvalue, rv.theta );

    //
    // Computing error bounds.
    //
    rv.log_error_bound = new double[ L * Q ];
    for( int i = 0; i < L * Q; i++ )
        rv.log_error_bound[ i ] = log(
            5 * L * ( N * K * log(
                ( double )N ) / log(
                2.0 ) + 20 * log(
                ( double )L * Q ) ) / log(
                2.0 ) ) + log(
            EPS ) - rv.theta * rv.step * i + L * rv.log_mgf;

    compute_theta_bounds(
        rv.log_error_bound, rv.step, L * Q, s, L * K, rv.log_min_error_bound,
        rv.log_max_error_bound, rv.theta );

    //
    // We estimate the p-value to be returned to be the geometic mean of the bounds.
    //
    rv.log_pvalue = ( rv.log_max_pvalue + rv.log_min_pvalue ) / 2;

    return rv;
}

//
// Computing the p-value for the entropy score using the csFFT algorithm.
// 
entropy_return_value
HS_csFFT(
    int N, int K, int L, double s, double* pu, int Q ) {

    //
    // rv = Value to be returned.
    //
    entropy_return_value rv;
    rv.Q = Q;

    //
    // Computing the pmf for a single column.
    //
    int* index;
    int size;
    double* log_pmf = HSlist(
        N, K, pu, Q, rv.step, index, size );

    //
    // Shifted L-fold convolution.
    //
    rv.theta = NOT_COMPUTED;
    rv.log_pmfL = csFFT(
        log_pmf, L, s, Q, rv.step, rv.theta, rv.log_mgf );
    compute_theta_pvalues(
        rv.log_pmfL, rv.step, L * Q, s, L * K, rv.log_min_pvalue,
        rv.log_max_pvalue, rv.theta );

    //
    // The current version does not compute error bounds for csFFT.
    //
    rv.log_error_bound = 0;

    //
    // We estimate the p-value to be returned to be the geometic mean of the bounds.
    //
    rv.log_pvalue = ( rv.log_max_pvalue + rv.log_min_pvalue ) / 2;

    return rv;
}

//
// Computing the p-value for the entropy score using the refine-csFFT algorithm.
//    accuracy = Desired number of decimal places of accuracy for the p-value
// 
entropy_return_value
HS_refine_csFFT(
    int N, int K, int L, double s, double* pu, int Q, double accuracy ) {

    //
    // rv = Value to be returned.
    //
    entropy_return_value rv;
    rv.Q = Q;

    //
    // Computing the pmf for a single column.
    //
    int* index;
    int size;
    double* log_pmf = HSlist(
        N, K, pu, Q, rv.step, index, size );

    //
    // Shifted L-fold convolution.
    //
    rv.theta = NOT_COMPUTED;
    rv.log_pmfL = new double[ L * Q ];
    refine_csFFT(
        log_pmf, L, s, Q, rv.step, accuracy, rv.theta, rv.log_mgf,
        rv.log_max_pvalue, rv.log_min_pvalue );

    //
    // The current version does not compute error bounds for csFFT.
    //
    rv.log_error_bound = 0;

    //
    // We estimate the p-value to be returned to be the geometic mean of the bounds.
    //
    rv.log_pvalue = ( rv.log_max_pvalue + rv.log_min_pvalue ) / 2;

    return rv;
}
