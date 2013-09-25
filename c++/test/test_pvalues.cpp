/**
 * Copyright John Reid 2011, 2012
 *
 * \file
 * \brief Code to test various p-value algorithms.
 *
 */

#define BOOST_TEST_MODULE compare pvalues test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/shared_array.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp>

#include <pvalues/calculator.h>

#include <iostream>

using namespace pvalues;
using namespace std;


BOOST_AUTO_TEST_CASE( compare )
{
    const double b[] = { .2, .2, .3, .3 };
    const std::vector< unsigned > ns = boost::assign::list_of
        // ( 1 )( 2 )
        // ( 5 )( 10 )( 20 )
        ( 100 )( 20 )
        // ( 400 )( 200 )
        // ( 800 )
        // ( 3000 )( 5000 )( 7000 )( 10000 )
        ;

    // create our Bejerano algorithm
    llr_pvalue_calculator::shared_ptr bejerano_algo = create_bejerano_pvalue_calculator( b, boost::size( b ), *std::max_element( ns.begin(), ns.end() ) );

    BOOST_FOREACH( unsigned n, ns ) {
        const size_t Q = 2 * n;

        // create our Hirji shifted algorithm
        boost::timer t;
        llr_pvalue_calculator::shared_ptr hs_algo = create_shifted_hirji_pvalue_calculator( b, boost::size( b ), n, Q );
        cout << "Ours took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << n << " with Q=" << Q << "\n";

        // create FAST algorithm
        t.restart();
        llr_pvalue_calculator::shared_ptr fast_algo = create_fast_pvalue_calculator( b, boost::size( b ), n, Q );
        cout << "FAST took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << n << " with Q=" << Q << "\n";

        const double delta = calculate_delta( b, boost::size( b ), n, Q );
        size_t ours_better_count = 0;
        size_t fast_better_count = 0;
        for( size_t q = 0; Q != q; ++q ) {
            const double LLR = q * delta;
            const double our_pval = hs_algo->operator()( n, LLR );
            const double fast_pval = fast_algo->operator()( n, LLR );
            boost::test_tools::check_is_close_t check_is_close;
            if( ! check_is_close( our_pval, fast_pval, boost::test_tools::percent_tolerance_t< double >( 1e-13 ) ) ) {
                const double bejerano_pval = bejerano_algo->operator()( n, LLR );
                const double our_diff = bejerano_pval - our_pval;
                const double fast_diff = bejerano_pval - fast_pval;
                if( std::fabs( our_diff ) < std::fabs( fast_diff ) ) {
                    ++ours_better_count;
                } else {
                    ++fast_better_count;
                }
//                cout
//                    << n << " : "
//                    << q << " : "
//                    << LLR << " : "
//                    << bejerano_pval << " : "
//                    << our_diff << " : "
//                    << fast_diff << " : "
//                    << ( std::fabs( our_diff ) < std::fabs( fast_diff ) ? "ours better" : "FAST better" )
//                    << "\n"
//                    ;
            }
            // BOOST_CHECK_CLOSE( our_pval, fast_pval, 1e-11 );
        }
        cout << "Ours better " << ours_better_count << " times, FAST better " << fast_better_count << " times.\n";
        BOOST_CHECK_GE( ours_better_count, fast_better_count );
    }
}
