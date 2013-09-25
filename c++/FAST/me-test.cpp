#include "utils.h"
#include "motif_evaluator.h"

//
// Test program for testing the motif_evaluator package.
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

    motif_evaluator m;

    //
    // Initializing the motif_evaluator object
    // and setting the number of sequences to N.
    //
    m.init(
        N, K, L, L, pu, Q );
    m.setN(
        N );

    //
    // Computing the p-value.
    //
    cout << setprecision(
        PRECISION );
    double log_pvalue = m.get_pvalue_careful(
        0, s );
    cout << "p-Value = " << log_pvalue << " (" << exp(
        log_pvalue ) << ")" << endl;

    return 0;
}
