#include <math.h>
#include <iostream>
#include "HS.h"

using namespace std;

#include "llr_score.h"
#include "multinom_prob.h"
#include "list.h"
#include "basic_pruning.h"
#include "multi_pmf.h"
#include "HS-options.h"
#include "HS-algorithm.h"

//
// Function: Computes and returns the entire lattice pmf (log-values in an array of
//           size Q) for the multinomial llr score based on the shifted-hirji algorithm
//           using lists.
//
// Parameters: N = Number of objects 
//             K = Number of bins 
//             pu = Null multinomial distribution 
//             Q = Lattice size 
//             step = Step size for the lattice (set by the function)
//             index = Array of indices for non-zero pmf entries (set by the function)
//             size = Size of the index array (set by the function)
//
double*
HSlist(
    int N, int K, double *pu, int Q, double& step, int*& index, int& size ) {

    array_return_t result = HSalgo< llr_score_t, multinom_prob_t< double > ,
                    list_t< double > , basic_pruning_t< llr_score_t,
                                    multinom_prob_t< double > , multi_pmf_t<
                                                    list_t< double > > > ,
                    array_return_t > (
        N, K, 0, pu, Q, set_shifted_p, get_array_return_value );

    step = result.step;
    index = result.index;
    size = result.size;
    return result.log_pmf;
}
