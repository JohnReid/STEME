#ifndef _UTILS_H_
#define _UTILS_H_

#define LOGZERO -1e100 
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define DEBUG 0
#define PRECISION 20

#include <iomanip>
#include <iostream>
#include "extended_exponent.h"

using namespace std;

//
// Double machine precision
//
#define EPS 2.220446049250313e-16

//
// Rounds to nearest integer
//
#define NINT(x) ((int) ((x)+((x) > 0 ? 0.5 : -0.5)))

#define LOG2(X) (log((double) X)/log(2.0))
#define RND_UP_TO_POWER_2(Y) (01 << (int) ceil(LOG2(Y)))

const double pi = 3.14159265358979;

//
// Computes log(exp(log_a) + exp(log_b))
//
inline double
log_sum(
    double log_a, double log_b ) {

    return ( ( log_a > log_b )
        ? log_a + log(
            1 + exp(
                log_b - log_a ) ) : log_b + log(
            1 + exp(
                log_a - log_b ) ) );
}

//
// Computes log_sum for the entire array of size n
//
double
log_array_sum(
    double* log_x, int n );

//
// Returns data[min_index] where data[min_index] = min(data[0..(size-1)])
//  
double
min_array(
    double* data, int size, int& min_index );

//
// Sets out[0..(row_size-1)][0..(column_size-1)] to
// exp(data[0..(row_size-1)][0..(column_size-1)]).
//
void
exp_2Darray(
    double** data, int row_size, int column_size, double** out );
void
exp_2Darray(
    double** data, int row_size, int column_size, ee_t** out );

inline int
sign(
    double x ) {
    return ( x >= 0
        ? 1 : -1 );
}
inline double
log(
    int x ) {
    return log(
        ( double )x );
}

inline unsigned long long
mod(
    int a, int b, int c ) {
    return ( ( ( unsigned long long )a ) * b ) % c;
}

inline double
log_exp_minus_one(
    double log_x ) {
    return ( log_x > 700
        ? log_x : log_x > 0
            ? log(
                exp(
                    log_x ) - 1 ) : log(
                1 - exp(
                    log_x ) ) );
}

template< class t >
    t
    exp_minus_one(
        t x ) {
        return ( abs(
            x ) < 1
            ? sinh(
                x / 2.0 ) * exp(
                x / 2.0 ) * 2.0 : exp(
                x ) - 1.0 );
    }
//{ return exp(x)-1.0; }

//
// Computes the L1 norm of the vector x of size n
// 
double
one_norm(
    double* x, int n );

//
// Computes the L2 norm of the vector x of size n
// 
double
two_norm(
    double* x, int n );

//
// Reads the parameter list for a program that computes
// the p-value for the multinomial llr statistic
//
void
read_llr_parameters(
    int argc, char** argv, int& N, int& K, double& s, int& Q, double*& pu );

//
// Reads the parameter list for a program that computes
// the p-value for the entropy score
//
void
read_entropy_parameters(
    int argc, char** argv, int& N, int& K, int& L, double& s, int& Q,
    double*& pu );

//
// Computes the p-value bounds for the given score s
//
void
compute_pvalues(
    double* log_pmf, double step, int Q, double s, int K, double& log_min_pval,
    double& log_max_pval );

//
// Computes the p-value bounds for the given score s and theta
//
void
compute_theta_pvalues(
    double* log_pmf, double step, int Q, double s, int K, double& log_min_pval,
    double& log_max_pval, double theta );

//
// Computes the error bounds for the given score s and theta
//
void
compute_theta_bounds(
    double* log_pmf, double step, int Q, double s, int K, double& log_min_pval,
    double& log_max_pval, double theta );

#include "complex.h"

//
// Computes the p-value directly from the characterisitic function
//
void
compute_phi_pvalue(
    complex* phi, double step, int size, double start, int end,
    int uncertainity, double& log_min_pval, double& log_max_pval, double theta,
    double log_mgf );

#endif 

