#include "utils.h"
#include "entropy-combined.h"

//
// Test program for testing functions for computing the p-value
// of the entropy score.
//
int
main(
    int argc, char** argv ) {

    //
    // Reading input parameters from the command line
    //
    int N;
    int K;
    int L;
    double s;
    int Q;
    double* pu;
    read_entropy_parameters(
        argc, argv, N, K, L, s, Q, pu );

#ifdef HS_SFFT
    entropy_return_value result = HS_sFFT(N, K, L, s, pu, Q);
#endif

#ifdef HS_CSFFT
    entropy_return_value result = HS_csFFT(N, K, L, s, pu, Q);
#endif

#ifdef HS_REFINE_CSFFT
    entropy_return_value result = HS_refine_csFFT(N, K, L, s, pu, Q, 1);
#endif

    cout << setprecision(
        PRECISION );
    cout << "p-Value = " << result.log_pvalue << " (" << exp(
        result.log_pvalue ) << ")" << endl;

#if DEBUG
    cout << "Max = " << result.log_max_pvalue << " (" << exp(result.log_max_pvalue) << ")" << endl;
    cout << "Min = " << result.log_min_pvalue << " (" << exp(result.log_min_pvalue) << ")" << endl;
#endif

    return 0;
}
