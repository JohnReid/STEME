/**
 * Copyright John Reid 2011
 *
 * \file Code to test PVAL p-value calculations
 */

#include "pvalue_test_defs.h"

#include <boost/assign/list_of.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>

#include <iostream>

using namespace boost;
using namespace std;


p_value_test_case
test_case( int N ) {
    p_value_test_case test_case( N );
    test_case.ICs.push_back( 1.146875 / 2. );
    return test_case;
}

int
main( int argc, char * argv[] ) {
	int pval_algorithm = 3;
	const double pu[] = { .3, .3, .2, .2 };
	const bool use_qfast = true;

	pvalue_calculator::ptr pval_calc = create_PVAL_pvalue_calculator( pu, pval_algorithm, use_qfast );

	p_value_test_case_vec test_cases = assign::list_of
        ( test_case( 128   ) )
        ( test_case( 1280  ) )
        ( test_case( 12800 ) )
        ( test_case( 12800 ) )
        ( test_case( 128   ) )
        ( test_case( 1280  ) )
        ( test_case( 12800 ) )
        ;

	BOOST_FOREACH( const p_value_test_case & tc, test_cases ) {
	    timer t;
	    const double log_pvalue = pval_calc->operator()( tc );
	    cout << "Took " << t.elapsed() << " seconds to calculate log p-value = " << log_pvalue << " for test case: " << tc << "\n";
	}

	return 0;
}
