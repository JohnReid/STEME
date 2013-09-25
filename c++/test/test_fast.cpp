/**
 * Copyright John Reid 2011
 *
 * \file
 * \brief Code to test FAST implementation of Hirji algorithm with Nagarajan's shifts.
 *
 */

#define BOOST_TEST_MODULE FAST pvalues test
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

// FAST HS algorithm
#include <HS.h>

using namespace pvalues;
using namespace std;





void
check_FAST( size_t N, size_t Q, const double * pu, size_t K ) {
    double step;
    int * index;
    int size;

    // run FAST HSlist algorithm
    boost::timer t;
    double * fast_log_pmf = HSlist( N, K, const_cast< double * >( pu ), Q, step, index, size );
    cout << "FAST took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << N << " with Q=" << Q << "\n";

    // check all indices are within Q
    for( int s = 0; size != s; ++s ) {
        // cout << index[ s ] << " : " << log_pmf[ index[ s ] ] << "\n";
        BOOST_CHECK( index[ s ] < int( Q ) );
    }

    // create FAST cumulative mass function
    boost::shared_array< double > fast_cmf = fast_cumulative( fast_log_pmf, Q, index, size, false );
    BOOST_CHECK_SMALL( fast_cmf[ 0 ], 6e-12 ); // check p.m.f. adds to 1.

    // normalise c.m.f.
    boost::iterator_range< double * > fast_cmf_range = boost::make_iterator_range( fast_cmf.get(), fast_cmf.get() + Q );
    normalise_log_cumulative( fast_cmf_range );
}




BOOST_AUTO_TEST_CASE( fast_1 )
{
    const double b[] = { .2, .2, .3, .3 };
    const std::vector< unsigned > ns = boost::assign::list_of( 1 )( 2 )( 5 )( 10 )( 20 );

    BOOST_FOREACH( unsigned n, ns ) {
        const size_t Q = n * 5 / 3 + 1;
        check_FAST( n, Q, b, boost::size( b ) );
    }
}


#ifdef NDEBUG

BOOST_AUTO_TEST_CASE( FAST )
{
    const double b[] = { .2, .2, .3, .3 };
    //const double b[] = { .3, .3, .2, .2 };
    const unsigned n = 1000;
    const unsigned Q = 2000;

    double step;
    int * index;
    int size;

    boost::timer t;
    double* log_pmf = HSlist( n, boost::size( b ), const_cast< double * >( b ), Q, step, index, size );
    cout << "FAST took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << n << " with Q=" << Q << "\n";
    double sum = 0.;
    for( unsigned q = 0; Q != q; ++q ) {
        //cout << q << " : " << log_pmf[ q ] << "\n";
        sum += std::exp( log_pmf[ q ] );
    }
    cout << "sum = " << sum << "\n";
}

#endif //NDEBUG

