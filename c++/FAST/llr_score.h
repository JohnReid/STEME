#ifndef _LLR_SCORE_H_
#define _LLR_SCORE_H_

#include "utils.h"

//
// Class for encapsulating the data and the methods for accessing them
// for the latticed log-likelihood ratio score (for use in Hirji's
// algorithm and bagFFT.) If a different score, say the chi-squared
// score needs to be used, then we can do so by creating a class which
// has the same public interface as in this class.
//
class llr_score_t {

private:

    //
    // N = Number of objects
    // K = Number of bins
    // pu = Background probabilities
    // Q = Lattice size
    // Ip[k][n] = Latticed llr score contribution for the kth bin having
    //            n objects in it
    // I[k][n] = llr score contribution for the kth bin having
    //           n objects in it
    //
    int N;
    int K;
    double* pu;
    int Q;
    int** Ip;
    double** I;

    //
    // step = Step size for the lattice
    // max_ent = Maximum lattice score
    // min_t_ent = Minimum partial lattice score
    //
    // Note: All lattice scores, even partial ones should lie
    //       between min_t_ent and max_ent.
    //
    double step;
    double max_ent;
    double min_ent;
    double min_t_ent;

public:

    //
    // Initializes all the data members
    // Q = Lattice size
    //
    llr_score_t(
        int given_N, int given_K, double* given_pu, int given_Q );

    inline int
    get_Ip(
        int k, int n ) {
        return Ip[ k ][ n ];
    }
    inline double
    get_I(
        int k, int n ) {
        return I[ k ][ n ];
    }

    //
    // Note: The caller should not modify the "copy".
    //
    inline int**
    get_Ip_copy() {
        return Ip;
    }

    inline int
    get_N() {
        return N;
    }
    inline int
    get_K() {
        return K;
    }
    inline int
    get_Q() {
        return Q;
    }
    inline double
    get_step() {
        return step;
    }
    inline double
    get_min_t_ent() {
        return min_t_ent;
    }
    inline double
    get_max_ent() {
        return max_ent;
    }
    inline double
    get_min_ent() {
        return min_ent;
    }

    ~llr_score_t();
};

#endif

