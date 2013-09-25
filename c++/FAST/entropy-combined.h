#ifndef _ENTROPY_COMBINED_H_
#define _ENTROPY_COMBINED_H_

//
// Type of return value for functions computing
// the p-value of the entropy score.
//
struct entropy_return_value {

    //
    // An estimate for the p-value;
    //
    double log_pvalue;

    //
    // Upper and lower bounds for the p-value and
    // the corresponding bounds on the absolute error.
    //
    double log_max_pvalue;
    double log_min_pvalue;
    double log_max_error_bound;
    double log_min_error_bound;

    //
    // step = Step size for the lattice
    // theta = Optimal shift
    // log_mgf = Log of the mgf at theta
    //
    double step;
    double theta;
    double log_mgf;

    //
    // Lattice pmf and the corresponding absolute error bounds
    // (log-values in an array of size Q).
    //
    int Q;
    double* log_pmfL;
    double* log_error_bound;
};

//
// Computing the p-value for the entropy score using the sFFT algorithm.
// 
entropy_return_value
HS_sFFT(
    int N, int K, int L, double s, double* pu, int Q );

//
// Computing the p-value for the entropy score using the csFFT algorithm.
//
entropy_return_value
HS_csFFT(
    int N, int K, int L, double s, double* pu, int Q );

//
// Computing the p-value for the entropy score using the refine-csFFT algorithm.
//    accuracy = Desired number of decimal places of accuracy for the p-value
// 
entropy_return_value
HS_refine_csFFT(
    int N, int K, int L, double s, double* pu, int Q, double accuracy );

#endif
