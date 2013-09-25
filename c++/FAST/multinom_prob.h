#ifndef _MULTINOM_PROB_H_
#define _MULTINOM_PROB_H_

//
// Class for encapsulating the data and the methods for accessing them
// for the multinomial probability distribution (for use in Hirji's
// algorithm and bagFFT.) If a different distribution needs to be
// used, then we can do so by creating a class which has the same
// public interface as in this class.
//
// p_t = Data type for storing the probability contributions
//
template< class p_t >
    class multinom_prob_t {

    private:

        //
        // N = Number of objects
        // K = Number of bins
        // p[k][n] = multinomial probability contribution for the kth bin
        //           having n objects in it = pu[k]^n/n!
        // log_p[k][n] = log(p[k][n])
        // log_factor[n] = log(n!)
        //
        int N;
        int K;
        double** log_p;
        p_t** p;
        double* log_factor;
        double* log_sum_pu;

        //
        // Thetas for the shifts
        //
        double theta1;
        double theta2;

    public:

        //
        // Initializes all the data members
        //
        multinom_prob_t(
            int given_N, int given_K, double* pu );

        //
        // Applies the theta1 shift to log_p
        //
        void
        shift1(
            double theta, int** Ip, double step );
        void
        norm1(
            double norm_factor );

        //
        // Applies the theta2 shift to log_p
        //
        void
        shift2(
            double theta );
        void
        norm2(
            int k, double norm_factor );

        void
        set_p() {
            exp_2Darray(
                log_p, K, N + 1, p );
        }
        inline double
        get_log_p(
            int k, int n ) {
            return log_p[ k ][ n ];
        }
        inline double**
        get_log_p_copy() {
            return log_p;
        }
        inline p_t
        get_p(
            int k, int n ) {
            return p[ k ][ n ];
        }
        inline p_t**
        get_p_copy() {
            return p;
        }
        inline int
        get_N() {
            return N;
        }
        inline int
        get_K() {
            return K;
        }
        inline double
        get_theta1() {
            return theta1;
        }
        inline double
        get_theta2() {
            return theta2;
        }

        //
        // Returns log of the missing constant factor in the probability
        // i.e. log(N!)
        //
        inline double
        log_const_factor() {
            return log_factor[ N ];
        }

        inline double
        log_factorial(
            int i ) {
            return log_factor[ i ];
        }

        inline double
        log_const_rem_factor(
            int k, int n ) {
            return log_factor[ N ] - log_factor[ N - n ] + ( N - n )
                            * log_sum_pu[ k ];
        }

        ~multinom_prob_t();
    };

#include "multinom_prob.cpp"

#endif
