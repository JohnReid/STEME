#ifndef _HS_ALGORITHM_H_
#define _HS_ALGORITHM_H_

#include <math.h>
#include "multi_pmf.h"

//
// Function: Implements hirji's DP algorithm.
//
// Parameters: N = Number of objects 
//             K = Number of bins 
//             s = Score
//             pu = Null distribution 
//             Q = Lattice size 
//             set_p = Function to set probabilities
//                     for computation in either extended-exponent
//                     or shifted mode (see HS-options.h)
//             get_return_value = Function to extract the return
//                                value in an appropriate format:
//                                array or single (see HS-options.h)
//
// score_t = Score type (see llr_score.h)
// prob_t = Probability type (see multinom_prob.h)
// pmf_t = Type for data structure to store the pmf (see list.h and array.h)
// prune_t = Pruning type (see basic_pruning.h and adv_pruning.h)
// return_t = Return type (see options.h)
//
template< class score_t, class prob_t, class pmf_t, class prune_t,
                class return_t >
    inline return_t
    HSalgo(
        int N, int K, double s, double *pu, int Q, void
        (*set_p)(
            prob_t& prob, score_t& score ), return_t
        (*get_return_value)(
            multi_pmf_t< pmf_t >& multi_pmf, int N, double step, int offset,
            prob_t& prob, prune_t& pruning ) ) {

        score_t score(
            N, K, pu, Q );
        prob_t prob(
            N, K, pu );

        set_p(
            prob, score );

        //
        // offset = Offset to latticed llr score to ensure that all array indices are positive
        // max_range = Maximum index for latticed llr score (after the offset has been applied)
        //
        int offset = ( int )( score.get_min_t_ent() / score.get_step() - K );
        int max_range = Q - offset;

        multi_pmf_t< pmf_t > prev_multi_pmf(
            N + 1, max_range );
        prune_t pruning(
            &score, &prob, &prev_multi_pmf, offset, s );

        for( int n = 0; n < N + 1; n++ )
            prev_multi_pmf.pmf[ n ].set(
                score.get_Ip(
                    0, n ) - offset, prob.get_p(
                    0, n ) );

        //
        // kth stage of recursion
        //
        for( int k = 1; k < K; k++ ){

            multi_pmf_t< pmf_t > curr_multi_pmf(
                N + 1, max_range );

            //
            // n is the number of objects in the bins.
            //
            for( int n = N; n >= 0; n-- ){
                //
                // i is the number of objects in the kth bin.
                //
                for( int i = 0; i <= n; i++ )
                    for( int It = pruning.start(
                        k - 1, n, i ); pruning.not_last(
                        k - 1, n, i ); It = pruning.next(
                        k - 1, n, i ) )
                        if( It + score.get_Ip(
                            k, i ) >= 0 && It >= 0 )
                            curr_multi_pmf.pmf[ n ].set(
                                It + score.get_Ip(
                                    k, i ), curr_multi_pmf.pmf[ n ][ It
                                                + score.get_Ip(
                                                    k, i ) ] + prob.get_p(
                                    k, i ) * prev_multi_pmf.pmf[ n - i ][ It ] );
            }

            prev_multi_pmf.destructive_copy(
                curr_multi_pmf );
        }

#if DEBUG
        cout << "Qavg " << pruning.get_Qavg() << endl;
#endif

        return get_return_value(
            prev_multi_pmf, N, score.get_step(), offset, prob, pruning );
    }

#endif
