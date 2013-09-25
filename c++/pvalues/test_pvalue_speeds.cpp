/**
 * Copyright John Reid 2011, 2013
 *
 * \file Code to test speed of various p-value calculation implementations
 */


#include <boost/timer.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>

#include "pvalue_test_defs.h"

using namespace boost;
using namespace std;
namespace po = boost::program_options;

void
apply_pvalue_calculator( const p_value_test_case_vec & arguments, pvalue_calculator::ptr calc, vector< double > & results ) {
	boost::function< double ( const p_value_test_case & ) > f = boost::ref( *calc );
	results.clear();
	results.reserve( arguments.size() );
	copy(
		make_transform_iterator( arguments.begin(), f ),
		make_transform_iterator( arguments.end(), f ),
		back_insert_iterator< vector< double > >( results )
	);
}

int
main( int argc, char * argv[] ) {

	size_t max_test_cases = 0;
    size_t max_N = 0;
    size_t min_N = 0;
#ifdef STEME_HAS_MEME_SOURCE
    bool do_meme = false;
#endif // STEME_HAS_MEME_SOURCE
	bool do_fast = false;
    bool do_pval = false;
    bool do_hs = false;
	bool use_qfast = false;
	int pval_algorithm = 3;
    size_t lattice_size = 0;

	//
	// Declare the supported options.
	//
	po::options_description generic_opts( "Allowed options" );
	generic_opts.add_options()
	    ( "help,h", "produce help message" )
	    ( "qfast", po::bool_switch( &use_qfast ), "use the MEME QFAST algorithm to combine column p-values" )
	    ( "pval", po::bool_switch( &do_pval ), "run PVAL on the p-value test cases" )
	    ( "fast", po::bool_switch( &do_fast ), "run FAST on the p-value test cases" )
#ifdef STEME_HAS_MEME_SOURCE
        ( "meme", po::bool_switch( &do_meme ), "run MEME on the p-value test cases" )
#endif // STEME_HAS_MEME_SOURCE
        ( "hs", po::bool_switch( &do_hs ), "run HS on the p-value test cases" )
	    ( "max-test-cases", po::value< size_t >( &max_test_cases ), "set limit on number of p-values to test" )
        ( "min-N", po::value< size_t >( &min_N ), "ignore p-values with low N" )
        ( "max-N", po::value< size_t >( &max_N ), "ignore p-values with high N" )
        ( "Q", po::value< size_t >( &lattice_size ), "lattice size" )
	    ( "pval-chi2", "Use the chi^2 approximate PVAL method" )
	    ( "pval-exact", "Use the exact PVAL method" )
	    ( "pval-bound", "Use the branch-and-bound PVAL method" )
	    ( "pval-convex", "Use the convex PVAL method" )
	;

	po::options_description hidden_opts( "Hidden options" );
	hidden_opts.add_options()
	    ( "input-file", po::value< string >(), "input file" )
	    ;

	po::options_description cmdline_opts;
	cmdline_opts.add( generic_opts ).add( hidden_opts );

	po::options_description visible_opts;
	visible_opts.add( generic_opts );

    po::positional_options_description p;
    p.add( "input-file", -1 );

	//
	// Parse the options
	//
    po::variables_map vm;
    store( po::command_line_parser( argc, argv ).options( cmdline_opts ).positional( p ).run(), vm );
    notify( vm );

	if( vm.count( "help" ) ) {
	    cout << visible_opts << "\n";
	    return 1;
	}

	if( vm.count( "max-test-cases" ) ) cout << "Only testing the first " << max_test_cases << " p-value test cases.\n";
    if( vm.count( "max-N" ) ) cout << "Ignoring p-value test cases with N higher than " << max_N << ".\n";
    if( vm.count( "min-N" ) ) cout << "Ignoring p-value test cases with N lower than " << min_N << ".\n";
	if( use_qfast ) cout << "Using QFAST to combine column p-values.\n";
	else cout << "Not using QFAST to combine column p-values.\n";

	if( vm.count( "pval-chi2" ) ) pval_algorithm = 0;
	if( vm.count( "pval-exact" ) ) pval_algorithm = 1;
	if( vm.count( "pval-bound" ) ) pval_algorithm = 2;
	if( vm.count( "pval-convex" ) ) pval_algorithm = 3;


	if( ! vm.count( "input-file" ) ) {
		cout << visible_opts << "\n";
		return -1;
	}

	// if not asked to do any of the methods, do them all
	if( ! (
#ifdef STEME_HAS_MEME_SOURCE
        do_meme ||
#endif // STEME_HAS_MEME_SOURCE
        do_fast || do_pval || do_hs )
    ) {
#ifdef STEME_HAS_MEME_SOURCE
        do_meme =
#endif // STEME_HAS_MEME_SOURCE
        do_fast = do_pval = do_hs = true;
	}


	//
	// Load p-value arguments
	//
	p_value_test_case_vec arguments;
	ifstream input( vm[ "input-file" ].as< string >().c_str() );
	read_p_value_test_cases( input, arguments, min_N, max_N, max_test_cases );
	cout << "Read " << arguments.size() << " p-value test cases.\n";
	if( max_test_cases && arguments.size() > max_test_cases ) {
		arguments.resize( max_test_cases );
	}
	cout << "Testing " << arguments.size() << " p-value test cases.\n";

	//
	// Set background distribution
	//
	double pu[] = { 0.3, 0.2, 0.2, 0.3 };
	cout
		<< "Using background distribution:"
		<< " p(A)=" << pu[0]
		<< " p(C)=" << pu[1]
		<< " p(G)=" << pu[2]
		<< " p(T)=" << pu[3]
		<< "\n"
		;


	//
	// Find smallest/biggest arguments
	//
	size_t minN = numeric_limits< size_t >::max();
	size_t maxN = 0;
	size_t minW = numeric_limits< size_t >::max();
	size_t maxW = 0;
	double minIC = numeric_limits< double >::max();
	double maxIC = 0.;
	BOOST_FOREACH( p_value_test_case args, arguments ) {
		minN = std::min( args.N, minN );
		maxN = std::max( args.N, maxN );
		const size_t W = args.W();
		minW = std::min( W, minW );
		maxW = std::max( W, maxW );
		const double IC = args.IC();
		minIC = std::min( IC, minIC );
		maxIC = std::max( IC, maxIC );
	}
	cout << "N range is [" << minN << "," << maxN << "]\n";
	cout << "W range is [" << minW << "," << maxW << "]\n";
	cout << "IC range is [" << minIC << "," << maxIC << "]\n";

    //
    // Check lattice size
    //
    if( ! lattice_size ) {
        lattice_size = 2 * maxN;
    }
    cout << "Lattice size = " << lattice_size << "\n";

	timer t;
	typedef vector< double > pvalue_vec;
	typedef boost::shared_ptr< pvalue_vec > pvalue_vec_ptr;

	//
	// Time MEME initialisation
	//
//	const vector< int > max_sites = boost::assign::list_of(32)(64)(128)(256)(512)(1024)(2048);
//	BOOST_FOREACH( int ms, max_sites ) {
//		init_MEME_llr_pv_tables( 2, ms, pu );
//	}


    //
    // Time HS algorithm
    //
    pvalue_vec_ptr hs_pvalues;
    if( do_hs ) {
        cout << "***************** Testing HS p-value calculations ************************\n";
        hs_pvalues.reset( new vector< double > );
        t.restart();
        apply_pvalue_calculator( arguments, create_HS_pvalue_calculator( maxN, pu, lattice_size, use_qfast ), *hs_pvalues );
        cout << "HS took " << t.elapsed() << " seconds to generate " << hs_pvalues->size() << " p-values\n";
        cout << "**************************************************************************\n";
    }


	//
	// Time pval algorithm
	//
	pvalue_vec_ptr pval_pvalues;
	if( do_pval ) {
		cout << "***************** Testing PVAL p-value calculations **********************\n";
		pval_pvalues.reset( new vector< double > );
		t.restart();
		apply_pvalue_calculator( arguments, create_PVAL_pvalue_calculator( pu, pval_algorithm, use_qfast ), *pval_pvalues );
		cout << "PVAL took " << t.elapsed() << " seconds to generate " << pval_pvalues->size() << " p-values\n";
		cout << "**************************************************************************\n";
	}


	//
	// Time FAST algorithm
	//
	pvalue_vec_ptr fast_pvalues;
	if( do_fast ) {
		cout << "***************** Testing FAST p-value calculations **********************\n";
		fast_pvalues.reset( new vector< double > );
		t.restart();
		apply_pvalue_calculator( arguments, create_FAST_pvalue_calculator( pu, lattice_size, use_qfast ), *fast_pvalues );
		cout << "FAST took " << t.elapsed() << " seconds to generate " << fast_pvalues->size() << " p-values\n";
		cout << "**************************************************************************\n";
	}


#ifdef STEME_HAS_MEME_SOURCE
	//
	// Time MEME algorithm
	//
	pvalue_vec_ptr meme_pvalues;
	if( do_meme ) {
		cout << "***************** Testing MEME p-value calculations **********************\n";
		meme_pvalues.reset( new vector< double > );
		t.restart();
		apply_pvalue_calculator( arguments, create_MEME_pvalue_calculator( minN, maxN, pu, lattice_size, use_qfast ), *meme_pvalues );
		cout << "MEME took " << t.elapsed() << " seconds to generate " << meme_pvalues->size() << " p-values\n";
		cout << "**************************************************************************\n";
	}
#endif // STEME_HAS_MEME_SOURCE


	//
	// Write p-values to file
	//
	const char * output_filename = "p-values.out";
	cout << "Writing p-values to " << output_filename << "\n";
	ofstream output( output_filename );
	for( size_t i = 0; arguments.size() != i; ++i ) {
		output << arguments[ i ];
		if( pval_pvalues ) output << "; PVAL=" << ( *pval_pvalues )[ i ];
		if( fast_pvalues ) output << "; FAST=" << ( *fast_pvalues )[ i ];
#ifdef STEME_HAS_MEME_SOURCE
        if( meme_pvalues ) output << "; MEME=" << ( *meme_pvalues )[ i ];
#endif // STEME_HAS_MEME_SOURCE
		output << "\n";
	}


	//
	// Write output in format that pval will read
	//
	const char * pval_filename = "pval.in";
	cout << "Writing pval input file: " << pval_filename << "\n";
	ofstream pval_input( pval_filename );
	for( size_t i = 0; arguments.size() != i; ++i ) {
		const p_value_test_case & tc = arguments[i];
		for( size_t w = 0; tc.W() != w; ++w ) {
			pval_input
				<< "4 " << pu[0] << " " << pu[1] << " " << pu[2] << " " << pu[3]
				<< " 0 " << tc.N << " " << tc.N * 2 * tc.ICs[w] << " "
				<< pval_algorithm << " 1 2 3 4\n";
		}
	}
}



