#ifndef _BASIC_PRUNING_H_
#define _BASIC_PRUNING_H_

//
// Class for doing basic pruning of a multi_pmf_t (see multi_pmf.h)
// in Hirji's algorithm. The "pruning" based on this class
// can be used to generate the entire non-zero pmf. 
// For more aggressive pruning for computing a single p-value
// the class adv_pruning_t in adv_pruning.h (not included in this package) 
// should be used.
//
// score_t = The type of the score for which the pruning is being
//           computed 
// prob_t = The type of the probability distribution for which the pruning is being
//           computed 
// multi_pmf_t = The instantiated type of multi_pmf_t which
//               is being pruned
//
template< class score_t, class prob_t, class multi_pmf_t >
    class basic_pruning_t {

    private:

        //
        // Score and multi_pmf for which the pruning is being carried out
        //
        score_t* score;
        multi_pmf_t* multi_pmf;

        //
        // min_firstK_I[k][n] = Minimum possible partial lattice score
        //                      with n objects in the first k bins
        // max_firstK_I[k][n] = Maximum possible partial lattice score
        //                      with n objects in the first k bins
        //
        int** min_firstK_I;
        int** max_firstK_I;

        //
        // Qavg is the average number of lattice entries on average
        // that are processed. Typically this should be much lesser
        // than the original lattice size Q. kn2 is set to K*N*(N-1)/2
        // and is used to compute Qavg;
        //
        double Qavg;
        double kn2;

    public:

        //
        // s = Score for which pruning is carried out
        //
        basic_pruning_t(
            score_t* given_score, prob_t* given_prob,
            multi_pmf_t* given_multi_pmf, int offset, double s );

        //
        // Functions for iterating over the pmf.
        // The appropriate pmf functions with the computed bounds are
        // called.
        //
        inline int
        start(
            int k, int n, int i ) {
            return multi_pmf->pmf[ n - i ].start(
                min_firstK_I[ k ][ n - i ] );
        }
        inline int
        not_last(
            int k, int n, int i ) {
            return multi_pmf->pmf[ n - i ].not_last(
                max_firstK_I[ k ][ n - i ] );
        }
        inline int
        next(
            int k, int n, int i ) {
            Qavg += kn2;
            return multi_pmf->pmf[ n - i ].next(
                max_firstK_I[ k ][ n - i ] );
        }

        inline double
        get_Qavg() {
            return Qavg;
        }

        ~basic_pruning_t() {

            for( int k = 0; k < score->get_K(); k++ ){

                delete[] ( min_firstK_I[ k ] );
                delete[] ( max_firstK_I[ k ] );
            }

            delete[] ( min_firstK_I );
            delete[] ( max_firstK_I );
        }
    };

#include "basic_pruning.cpp"

#endif
