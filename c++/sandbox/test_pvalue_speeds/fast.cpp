/**
 * Copyright John Reid 2011
 *
 * \file Code to test speed of FAST p-value calculation implementation
 */


// FAST headers
#include "utils.h"
#include "motif_evaluator.h"

#include "pvalue_test_defs.h"

#include <boost/assert.hpp>
#include <boost/timer.hpp>


/**
 * A wrapper around Nagarajan's pvalue calculations
 */
struct fast_pvalue_calculator : pvalue_calculator, boost::noncopyable {
	typedef boost::shared_ptr< motif_evaluator > motif_evaluator_ptr;
	typedef std::map< int, motif_evaluator_ptr > motif_evaluator_map;

	motif_evaluator_map evaluators;
	int minL;
	int maxL;
	int Q;
	double * pu;

	fast_pvalue_calculator( int minL, int maxL, double * pu, int Q = 16384 )
	: minL( minL )
	, maxL( maxL )
	, Q( Q )
	, pu( pu )
	{ }

	double
	operator()( const p_value_test_case & args ) {
		// check everything is as it should be
		BOOST_ASSERT( int( args.W() ) >= minL );
		BOOST_ASSERT( int( args.W() ) <= maxL );

		return get_evaluator( args.N )->get_pvalue_careful( args.W() - minL, args.N * args.IC() );
	}


protected:
	motif_evaluator *
	get_evaluator( int N ) {
		motif_evaluator_map::iterator i = evaluators.find( N );
		if( evaluators.end() == i ) {
			cout << "Creating new evaluator for N=" << N << " ... ";
			bool inserted;
			boost::timer timer;
			boost::tie( i, inserted ) = evaluators.insert( motif_evaluator_map::value_type( N, motif_evaluator_ptr( new motif_evaluator ) ) );
			i->second->init( N, 4, minL, maxL, pu, Q );
			i->second->setN( N );
			BOOST_ASSERT( inserted );
			cout << "took " << timer.elapsed() << " seconds\n";
		}
		return i->second.get();
	}
};


pvalue_calculator::ptr
create_FAST_pvalue_calculator(
	int minL,
	int maxL,
	double * pu,
	int Q
) {
	return pvalue_calculator::ptr( new fast_pvalue_calculator( minL, maxL, pu, Q ) );
}







