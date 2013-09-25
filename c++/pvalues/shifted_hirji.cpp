/**
 * Copyright John Reid 2011, 2013
 *
 * \file Code to implement p-value calculation using shifted Hirji algorithm
 */


#include "pvalue_test_defs.h"

#include <pvalues/calculator.h>
#include <pvalues/hirji_shifted.h>

#include <iostream>

using namespace std;


namespace pvalues {

/**
 * A wrapper around Hirji shifted pvalue calculations
 */
struct hs_pvalue_calculator : pvalue_calculator {
	bool use_qfast;
    hirji_shifted<> alg;
    hirji_shifted<>::cumulative_distributions log_cmfs;

    hs_pvalue_calculator(
        int maxN,
        const double * pu,
        bool use_qfast,
        size_t Q
    )
    : use_qfast( use_qfast )
    , alg( boost::make_iterator_range( pu, pu + 4 ), maxN, Q )
    , log_cmfs( alg )
    {
        if( use_qfast ) cout << "Using QFAST to combine column p-values.\n";
        else cout << "Not using QFAST to combine column p-values.\n";
    }

    virtual ~hs_pvalue_calculator() { }

    double operator()( const p_value_test_case & args ) {
        // calculate the p-value
        double log_pop = 0.;
        BOOST_FOREACH( double IC, args.ICs ) {
            log_pop += log_cmfs.get_log_p_value( IC * args.N, args.N );
        }
        return use_qfast ? log_qfast( args.W(), log_pop ) : log_pop;
    }
};




/**
 * A p-value calculator that uses the shifted Hirji algorithm.
 */
struct shifted_hirji_pvalue_calculator
: llr_pvalue_calculator
, boost::noncopyable
{
    hirji_shifted< double > algorithm;
    boost::shared_array< double > log_cmf;

    shifted_hirji_pvalue_calculator( const double * bg, size_t K, size_t max_N, size_t Q )
    : algorithm(
        boost::make_iterator_range( bg, bg + K ),
        max_N,
        Q
    ) {
        // create cumulative density function
        hirji_shifted< double >::slice_t log_pmf = algorithm.get_log_pmf( max_N );
        log_cmf = calculate_log_cumulative( log_pmf );
    }

    /// Calculate the p-value of a motif column's LLR statistic
    double
    operator()( size_t N, double LLR ) {
        if( N > algorithm.n ) {
            throw std::logic_error( "N too large." );
        }
        return log_cmf[ get_LLR_index( algorithm.delta, LLR ) ];
    }

    /// \return The largest N, this calculator can handle.
    size_t
    get_max_N() {
        return algorithm.n;
    }
};



/**
 * Create a p-value calculator based on the shifted Hirji algorithm.
 */
boost::shared_ptr< llr_pvalue_calculator >
create_shifted_hirji_pvalue_calculator( const double * bg, size_t K, size_t N, size_t Q ) {
    return boost::shared_ptr< llr_pvalue_calculator >( new shifted_hirji_pvalue_calculator( bg, K, N, Q ) );
}



} // namespace pvalues




pvalue_calculator::ptr
create_HS_pvalue_calculator(
    int maxN,
    const double * pu,
    size_t Q,
    bool use_qfast
) {
    return pvalue_calculator::ptr( new pvalues::hs_pvalue_calculator( maxN, pu, use_qfast, Q ) );
}


