#include "utils.h"
#include "naive_enumeration.cpp"
#include "convolution.h"
#include "motif_evaluator.h"

//
// Initializes the members of the class.
//
void
motif_evaluator::init(
    int given_maxN, int given_K, int given_minL, int given_maxL,
    double* given_pu, int given_Q ) {

    maxN = given_maxN;
    K = given_K;
    minL = given_minL;
    maxL = given_maxL;
    pu = given_pu;
    Q = given_Q;
    log_pmf = multi_naive_enumeration(
        maxN, pu, Q, step );

    log_pvalues = new double*[ maxL - minL + 1 ];
    log_acc = new double*[ maxL - minL + 1 ];
    for( int L = minL; L <= maxL; L++ ){
        log_pvalues[ L - minL ] = new double[ L * Q ];
        log_acc[ L - minL ] = new double[ L * Q ];
    }
}

#include "minimize.h"

//
// Class for optimizing the value of theta (see thetaFn.h for another
// example).
//
class thetaFn {

private:

    static llr_score_t* score;
    static multinom_prob_t< ee_t >* prob;
    static int N;
    static double st;

public:

    static void
    init(
        int given_N, int K, double* pu, int Q, double given_st ) {

        N = given_N;
        st = given_st;
        score = new llr_score_t(
            N, K, pu, Q );
        prob = new multinom_prob_t< ee_t > (
            N, K, pu );
        set_ee_p(
            *prob, *score );
    }

    static double
    minimizationFn(
        double theta ) {

        ee_t sum;
        for( int i = 0; i <= N; i++ )
            for( int j = 0; j <= N - i; j++ )
                for( int k = 0; k <= N - ( i + j ); k++ ){
                    int l = N - ( i + j + k );
                    double I = score->get_I(
                        0, i ) + score->get_I(
                        1, j ) + score->get_I(
                        2, k ) + score->get_I(
                        3, l );
                    ee_t p = prob->get_p(
                        0, i ) * prob->get_p(
                        1, j ) * prob->get_p(
                        2, k ) * prob->get_p(
                        3, l );
                    sum += p * ee_t(
                        theta * I );
                }

        return sum.log_get() + prob->log_const_factor() - theta * st;
    }

    static void
    clear() {

        delete ( score );
        delete ( prob );
    }
};

llr_score_t* thetaFn::score;
multinom_prob_t< ee_t >* thetaFn::prob;
int thetaFn::N;
double thetaFn::st;

//
// If s is given, the p-values are precomputed
// with an appropriate theta computed for it.
// Otherwise the given theta is used.
//  
void
motif_evaluator::precompute(
    double theta, double s ) {

    double* curr_log_pmf = log_pmf[ N ];
    double log_mgf = LOGZERO;

    int
                    curr_Q =
                                    MIN(Q, int(Q/((double) maxN)*N));
    int period = RND_UP_TO_POWER_2(maxL*curr_Q);

    //
    // Computing theta if needed.
    //
    if( s != -1 ){

        thetaFn::init(
            N, K, pu, curr_Q, s );
        log_mgf = minimize(
            thetaFn::minimizationFn, theta );
        log_mgf += theta * s;
        thetaFn::clear();

#if DEBUG
        cout << "Precomputed with: " << theta << " " << s << endl;
#endif
    }

    double* pmf = new double[ period ];
    for( int It = curr_Q; It < period; It++ )
        pmf[ It ] = 0;

    if( log_mgf == LOGZERO )
        log_mgf = log_array_sum(
            pmf, curr_Q );

    for( int It = 0; It < curr_Q; It++ )
        pmf[ It ] = exp(
            curr_log_pmf[ It ] + theta * step * It - log_mgf );

    //
    // Simultaneously computing all the pmfs.
    //
#if NR
    double** pmfL = multi_self_conv_fftnr(pmf, minL, maxL, period);
#else
    double** pmfL = multi_self_conv_fftw(
        pmf, minL, maxL, period );
#endif

    delete[] ( pmf );

    //
    // Updating the p-values and error bounds.
    //
    for( int L = minL; L <= maxL; L++ ){

        double log_err_term = log(
            EPS * L * ( 20 * log(
                period ) + 20 ) ) + L * log_mgf;
        double log_err_sum = LOGZERO;
        double log_pvalue_sum = LOGZERO;

        for( int i = L * curr_Q - 1; i < L * Q; i++ ){

            log_pvalues[ L - minL ][ i ] = -( L * curr_Q - 1 ) * step;
            log_acc[ L - minL ][ i ] = log_acc_level + 1;
        }

        if( theta >= 0 ){

            //log_pvalue_sum = -(L*curr_Q-1)*step;
            for( int It = L * curr_Q - 2; It >= 0; It-- ){

                if( pmfL[ L - minL ][ It ] > 0 ){
                    log_pvalue_sum = log_sum(
                        log_pvalue_sum, log(
                            pmfL[ L - minL ][ It ] ) - theta * step * It + L
                                        * log_mgf );
                    log_err_sum = log_sum(
                        log_err_sum, -theta * step * It );
                }

                if( log_pvalue_sum - ( log_err_sum + log_err_term )
                                > log_acc[ L - minL ][ It ] ){
                    log_pvalues[ L - minL ][ It ] = log_pvalue_sum;
                    log_acc[ L - minL ][ It ] = log_pvalue_sum - ( log_err_sum
                                    + log_err_term );

                }
            }
        }else{

            for( int It = 0; It < L * curr_Q - 1; It++ ){

                if( pmfL[ L - minL ][ It ] > 0 ){
                    log_pvalue_sum = log_sum(
                        log_pvalue_sum, log(
                            pmfL[ L - minL ][ It ] ) - theta * step * It + L
                                        * log_mgf );
                    log_err_sum = log_sum(
                        log_err_sum, -theta * step * It );
                }

                if( log(
                    1 - exp(
                        log_pvalue_sum ) ) - ( log_err_sum + log_err_term )
                                > log_acc[ L - minL ][ It ] ){
                    log_pvalues[ L - minL ][ It ] = log(
                        1 - exp(
                            log_pvalue_sum ) );
                    log_acc[ L - minL ][ It ] = log_pvalues[ L - minL ][ It ]
                                    - ( log_err_sum + log_err_term );
                }
            }
        }

        delete[] ( pmfL[ L - minL ] );
    }

    delete[] ( pmfL );
}

void
motif_evaluator::setN(
    int given_N ) {

    N = given_N;

    for( int L = minL; L <= maxL; L++ )
        for( int It = 0; It < L * Q; It++ ){
            log_pvalues[ L - minL ][ It ] = 0;
            log_acc[ L - minL ][ It ] = -100;
        }

    precompute(
        1 );
}

motif_evaluator::~motif_evaluator() {

    for( int i = 0; i < maxN + 1; i++ )
        delete[] ( log_pmf[ i ] );

    for( int i = minL; i <= maxL; i++ )
        delete[] ( log_pvalues[ i - minL ] );

    delete[] ( log_pmf );
    delete[] ( log_pvalues );
}

