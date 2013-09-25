/**
 * Copyright John Reid 2011, 2013
 *
 * \file Bejerano p-value calculations
 */


#include <pvalues/bejerano.h>
#include <pvalues/calculator.h>

namespace pvalues {



/**
 * A p-value calculator that uses Bejerano algorithm.
 */
struct bejerano_pvalue_calculator
: llr_pvalue_calculator
, boost::noncopyable
{
	virtual ~bejerano_pvalue_calculator() {}

    log_arithmetic<> log_arith;
    const std::vector< double > b;

    bejerano_pvalue_calculator( const double * bg, size_t K, size_t max_N )
    : log_arith( max_N )
    , b( bg, bg + K )
    {
    }

    /// Calculate the p-value of a motif column's LLR statistic
    double
    operator()( size_t N, double LLR ) {
        if( N > log_arith.max_n() ) {
            throw std::logic_error( "Bejerano calculator: N too large." );
        }
        bejerano<> algorithm( log_arith, N, LLR, b );
        return algorithm();
    }

    /// \return The largest N, this calculator can handle.
    size_t
    get_max_N() {
        return log_arith.max_n();
    }
};

/**
 * Create a p-value calculator that uses the Bejerano algorithm.
 */
boost::shared_ptr< llr_pvalue_calculator >
create_bejerano_pvalue_calculator( const double * bg, size_t K, size_t max_N ) {
    return boost::shared_ptr< llr_pvalue_calculator >( new bejerano_pvalue_calculator( bg, K, max_N ) );
}





} // namespace pvalues
