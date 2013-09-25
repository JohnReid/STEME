/**
 * Copyright John Reid 2011
 *
 * \file Code to calculate p-values.
 */


#ifndef _PVALUES_JR_25AUG2011_CALCULATOR_H_
#define _PVALUES_JR_25AUG2011_CALCULATOR_H_

#include <pvalues/defs.h>

namespace pvalues {




/**
 * Calculates the p-value of a motif column's LLR statistic.
 */
struct llr_pvalue_calculator {

    typedef boost::shared_ptr< llr_pvalue_calculator > shared_ptr;

    /// Virtual destructor.
    virtual
    ~llr_pvalue_calculator();

    /// Calculate the p-value of a motif column's LLR statistic.
    virtual
    double
    operator()( size_t N, double LLR ) = 0;

    /// \return The largest N, this calculator can handle.
    virtual
    size_t
    get_max_N() = 0;

};



/**********************************************************************/
/*
    log_qfast

    Calculate the log p-value of the log of the
    product of uniform [0,1] random variables.

    Derived from MEME source code.

*/
/**********************************************************************/
inline
double
log_qfast(
    int n,               /* number of random variables in product */
    double logk          /* product of random variables */
) {
    int i;
    double term, phi;

    if( n == 0 ) {
        return 0; /* worst possible log p-value */
    }

    phi = term = 1;
    for( i = 1; i != n; i++ ) {
        term *= -logk / i;
        phi += term;
    }

    return logk + std::log( phi );
}                               /* qfast */



/**********************************************************************/
/*
    get_log_nalign

    Get an upper bound on the number of independent alignments
    of segments of length w.

    Same as get_log_nalign() in meme source when model is tcm and invcomp is true.

    Derived from MEME source code.

*/
/**********************************************************************/
inline
double
get_log_nalign(
    int N,                      /* number of occurrences */
    int num_possible_sites      /* number of possible sites */
)
{
    double log_nalign = 0.;        /* log number alignments */

    if( N > num_possible_sites ) {                       /* impossible N */

        log_nalign = std::numeric_limits< double >::max();

    } else { /* remove 1 site per site */

        for( int i = 0; i != N; ++i ) {
            log_nalign += std::log( (num_possible_sites - i) * 2 / (i + 1) );
        }

    }

    return log_nalign;
}                               /* double get_log_nalign */







/**
 * Calculate the delta used to index LLRs in a lattice.
 */
double
calculate_delta( const double * bg, size_t K, size_t N, size_t Q );



/**
 * Get the index of a LLR in a lattice.
 */
size_t
get_LLR_index( double delta, double LLR );


/**
 * Create a cumulative distribution from a FAST p.m.f.
 */
boost::shared_array< double >
fast_cumulative(
    double * log_pmf,
    int Q,
    int * index,
    int size,
    bool normalise = true
);


/**
 * Create a p-value calculator that uses the FAST code.
 */
llr_pvalue_calculator::shared_ptr
create_fast_pvalue_calculator( const double * bg, size_t K, size_t N, size_t Q );



/**
 * Create a p-value calculator based on the shifted Hirji algorithm.
 */
llr_pvalue_calculator::shared_ptr
create_shifted_hirji_pvalue_calculator( const double * bg, size_t K, size_t N, size_t Q );


/**
 * Create a p-value calculator that uses the Bejerano algorithm.
 */
boost::shared_ptr< llr_pvalue_calculator >
create_bejerano_pvalue_calculator( const double * bg, size_t K, size_t max_N );





} // namespace pvalues

#endif /* _PVALUES_JR_25AUG2011_CALCULATOR_H_ */
