/**
 * Copyright John Reid 2011, 2013
 *
 * \file Code to test speed of FAST p-value calculation implementation
 */


// FAST headers
#include <HS.h>

#include "pvalue_test_defs.h"

#include <pvalues/calculator.h>

#include <boost/assert.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/math/special_functions/round.hpp>

using namespace std;

/**
 * A wrapper around Nagarajan's pvalue calculations
 */
struct fast_pvalue_calculator : pvalue_calculator, boost::noncopyable {

	typedef boost::shared_array< double > HS_array;
	typedef std::map< int, HS_array > HS_array_map;

	HS_array_map pmfs;
    const bool use_qfast;
	const int Q;
	const double * pu;
	const double max_IC;

	fast_pvalue_calculator( const double * pu, int Q = 16384, bool use_qfast = false )
	: use_qfast( use_qfast )
	, Q( Q )
	, pu( pu )
	, max_IC( - std::log( *std::min_element( pu, pu + 4 ) ) )
	{
        if( use_qfast ) cout << "Using QFAST to combine column p-values.\n";
        else cout << "Not using QFAST to combine column p-values.\n";
	}

	virtual ~fast_pvalue_calculator() { }

	size_t
	index_for_IC( double IC ) {
	    return boost::math::round( IC / max_IC * Q );
	}

	double
	operator()( const p_value_test_case & args ) {
        double log_pop = 0.;
        const double * array = get_array( args.N );
        BOOST_FOREACH( double IC, args.ICs ) {
            log_pop += array[ index_for_IC( IC ) ];
        }
        return use_qfast ? pvalues::log_qfast( args.W(), log_pop ) : log_pop;
	}


protected:
	double *
	get_array( int N ) {
	    HS_array_map::iterator i = pmfs.find( N );
		if( pmfs.end() == i ) {
			cout << "Creating new evaluator for N=" << N << " ... ";
	        double step;
	        int * index;
	        int size;
			bool inserted;
			boost::timer timer;
			boost::tie( i, inserted ) = pmfs.insert(
			    HS_array_map::value_type(
			        N,
			        HS_array(
			            HSlist(
			                N,
			                4,
			                const_cast< double * >( pu ),
			                Q,
			                step,
			                index,
			                size
			            )
			        )
                )
			);
			BOOST_ASSERT( inserted );
			cout << "took " << timer.elapsed() << " seconds\n";
		}
		return i->second.get();
	}
};


pvalue_calculator::ptr
create_FAST_pvalue_calculator(
	const double * pu,
	int Q,
	bool use_qfast
) {
	return pvalue_calculator::ptr( new fast_pvalue_calculator( pu, Q, use_qfast ) );
}






namespace pvalues {




llr_pvalue_calculator::~llr_pvalue_calculator() { }



/**
 * Calculate a cumulative density function
 */
boost::shared_array< double >
fast_cumulative(
    double * log_pmf,
    int Q,
    int * index,
    int size,
    bool normalise
) {
    boost::shared_array< double > result( new double[ Q ] );
    std::sort( index, index + size );
    int s = size - 1;
    double log_cdf = std::log( 0. );
    for( int q = Q; 0 != q; ) {
        --q;
        if( s >= 0 && q == index[ s ] ) {
            log_cdf = log_add( log_cdf, log_pmf[ q ] );
            --s;
            BOOST_ASSERT( s >= -1 );
        }
        result[ q ] = log_cdf;
    }
    if( normalise ) {
        boost::iterator_range< double * > cmf_range = boost::make_iterator_range( result.get(), result.get() + Q );
        const double norm_adjustment = normalise_log_cumulative( cmf_range );
        cout << "FAST normalisation adjustment " << norm_adjustment << "\n";
    }
    return result;
}


/**
 * A p-value calculator that uses the FAST code.
 */
struct fast_pvalue_calculator
: llr_pvalue_calculator
, boost::noncopyable
{
    size_t Q;
    size_t N;
    double step;
    int * index;
    int size;
    double * log_pmf;
    double delta;
    boost::shared_array< double > log_cmf;

    fast_pvalue_calculator( const double * bg, size_t K, size_t N, size_t Q )
    : Q( Q )
    , N( N )
    , log_pmf( HSlist( N, K, const_cast< double * >( bg ), Q, step, index, size ) )
    , delta( calculate_delta( bg, K, N, Q ) )
    {
        // create FAST cumulative mass function
        log_cmf = fast_cumulative( log_pmf, Q, index, size );
    }

    ~fast_pvalue_calculator() {
        delete [] log_pmf;
        delete [] index;
    }

    /// Calculate the p-value of a motif column's LLR statistic
    double
    operator()( size_t _N, double LLR ) {
        if( _N != N ) {
            throw std::logic_error( "FAST p-value calculator only works on one N." );
        }
        return log_cmf[ get_LLR_index( delta, LLR ) ];
    }

    /// \return The largest N, this calculator can handle.
    size_t
    get_max_N() {
        return N;
    }
};


/**
 * Calculate the delta used to index LLRs in a lattice.
 */
double
calculate_delta( const double * bg, size_t K, size_t N, size_t Q ) {
    const double log_min_bg = std::log( *std::min_element( bg, bg + K ) );
    const double max_LLR = - ( N * log_min_bg );
    return max_LLR / ( Q - 1 );
}


/**
 * Get the index of a LLR in a lattice.
 */
size_t
get_LLR_index( double delta, double LLR ) {
    return boost::math::round( LLR / delta );
}




/**
 * Create a p-value calculator that uses the FAST code.
 */
boost::shared_ptr< llr_pvalue_calculator >
create_fast_pvalue_calculator( const double * bg, size_t K, size_t N, size_t Q ) {
    return boost::shared_ptr< llr_pvalue_calculator >( new fast_pvalue_calculator( bg, K, N, Q ) );
}




} // namespace pvalues



