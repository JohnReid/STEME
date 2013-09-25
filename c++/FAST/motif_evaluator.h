#ifndef _MOTIF_EVALUATOR_H_
#define _MOTIF_EVALUATOR_H_

#include <vector>
using namespace std;

//
// Sufficient accuracy level for
// reporting p-values. Currently
// set to 2 decimal places.
//
const double log_acc_level = 2 * log(
    10 );

//
// Class for implementing the memo-sFFT
// algorithm.
//
class motif_evaluator {

private:

    //
    // maxN = Largest value of N (the depth of the motif)
    //        for which p-values are memoized
    // N = Current value of N for which p-values are reported
    // K = Size of the alphabet
    // minL = Minimum value for L (the width of the motif)
    // maxL = Maximum value for L (the width of the motif)
    //
    int maxN;
    int N;
    int K;
    int minL;
    int maxL;

    //
    // pu = Background frequencies for the letters of the alphabet
    // Q = Lattice size
    // step = Step size dictated by the lattice
    //
    double* pu;
    int Q;
    double step;
    double** log_pmf;

    //
    // log_pvalues[L-minL][I] = log p-value for a latticed entropy score of
    //                          I for a motif of width L
    // log_acc[L-minL][I] = log(p-value/error-bound) for the corresponding
    //                      p-value
    //
    double** log_pvalues;
    double** log_acc;

    //
    // If s is given, the p-values are precomputed
    // with an appropriate theta computed for it.
    // Otherwise the given theta is used.
    //
    void
    precompute(
        double theta, double s = -1 );

public:

    motif_evaluator() {
    }

    //
    // Initializes the members of the class.
    //
    void
    init(
        int given_maxN, int given_K, int given_minL, int given_maxL,
        double* given_pu, int given_Q );

    void
    setN(
        int given_N );
    int
    getN() {
        return N;
    }

    //
    // Tests to see if the memoized p-value is
    // accurate before it is returned.
    //
    double
    get_pvalue_careful(
        int L_minus_min, double s ) {

        int I = int(
            s / step );

        if( log_acc[ L_minus_min ][ I ] > log_acc_level )
            return log_pvalues[ L_minus_min ][ I ];
        else{

            precompute(
                0, s / ( ( double )L_minus_min + minL ) );
            return log_pvalues[ L_minus_min ][ I ];
        }
    }

    //
    // Returns the memoized p-value without testing
    // its accuracy.
    //
    double
    get_pvalue(
        int L_minus_min, double s ) {

        return log_pvalues[ L_minus_min ][ int(
            s / step ) ];
    }

    ~motif_evaluator();
};

#endif
