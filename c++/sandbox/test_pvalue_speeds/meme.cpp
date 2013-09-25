/**
 * Copyright John Reid 2011
 *
 * \file Code to test speed of MEME p-value calculation implementation
 */


#include "pvalue_test_defs.h"

extern "C" {
#include <steme/meme.h>
}

#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/timer.hpp>

#include <vector>

using namespace std;


/**
 * A wrapper around MEME's pvalue calculations
 */
struct meme_pvalue_calculator : pvalue_calculator {
	double * pu;
	int range;
	bool use_qfast;
	meme_pvalue_calculator(
		int minN,
		int maxN,
		double * pu,
		bool use_qfast,
		int range
	)
	: pu( pu )
	, range( range )
	, use_qfast( use_qfast )
	{
		if( use_qfast ) cout << "Using QFAST to combine column p-values.\n";
		else cout << "Not using QFAST to combine column p-values.\n";
		// initialise MEME's tables
		init_MEME_llr_pv_tables( minN, maxN, pu );
	}

	double operator()( const p_value_test_case & args ) {
		// calculate the p-value
		double log_pop = 0.;
		BOOST_FOREACH( double IC, args.ICs ) {
			log_pop += get_llr_pv( IC * args.N, args.N, 1, range, 1.0, 4, pu );
		}
		return use_qfast ? log_qfast( args.W(), log_pop ) : log_pop;
	}
};

pvalue_calculator::ptr
create_MEME_pvalue_calculator(
	int minN,
	int maxN,
	double * pu,
	int range,
	bool use_qfast
) {
	return pvalue_calculator::ptr( new meme_pvalue_calculator( minN, maxN, pu, use_qfast, range ) );
}




/**
 * Initialise MEME LLR pv tables.
 */
void
init_MEME_llr_pv_tables(
  int min,				              /* minimum number of sites */
  int max,				              /* maximum number of sites */
  double * pu                         /* Background frequencies. */
) {
	//initialise the tables.
	cout << "Initialising MEME's LLR p-value tables for numbers of sites from " << min << " to " << max << " ... ";
	init_log();
	init_exp();
	reset_llr_pv_tables();
	boost::timer timer;
	init_llr_pv_tables( min, max, 4, pu, 0 );
	cout << "took " << timer.elapsed() << " seconds\n";
}







