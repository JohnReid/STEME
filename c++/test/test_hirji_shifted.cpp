/**
 * Copyright John Reid 2011
 *
 * \file
 * \brief Code to test Hirji algorithm with Nagarajan's shifts.
 *
 */

#define BOOST_TEST_MODULE hirji shifted
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/shared_array.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp>

#include <pvalues/hirji_shifted.h>
#include <pvalues/calculator.h>

#include <iostream>

using namespace pvalues;
using namespace std;




void
check_is_proper_dist( hirji_shifted<> & alg, double percent_tol = 1e-9 ) {
    for( size_t n1 = 1; alg.n >= n1; ++n1 ) {
        double sum = 0.;
        int first_j, last_j = 0;
        const int j_begin = alg.get_minimal_j();
        const int j_end = alg.get_Q();
        for( int j = j_begin; j_end != j; ++j ) {
            // cout << j << " : " << p << "\n";
            const double p = std::exp( alg.get_p( j, n1 ) );
            if( p ) {
                if( ! sum ) {
                    first_j = j;
                }
                last_j = j;
            }
            sum += p;
        }
        // cout << n1 << " : " << first_j << " : " << last_j << "\n";
        BOOST_CHECK_EQUAL( first_j, alg.j_begin( n1 ) );
        // BOOST_CHECK_EQUAL( last_j + 1, alg.j_end( n1 ) );
        // cout << "sum = " << sum << "\n";
        BOOST_CHECK_CLOSE( 1., sum, percent_tol );
    }
}



template< typename InIt, typename OutIt >
double
log_cumulative( InIt begin, InIt end, OutIt out )
{
    double curr = std::log( 0. );
    while( begin != end ){
        curr = log_add( curr, *begin );
        ++begin;
        *out = curr;
        ++out;
    }
    return curr;
}

template< typename Range >
double
make_cumulative( const Range & log_pmf, std::vector< double > & out ) {
    typedef typename boost::range_value< Range >::type value_t;
    out.clear();
    return log_cumulative( boost::rbegin( log_pmf ), boost::rend( log_pmf ), std::back_inserter( out ) );
}


void
check_for_smaller_n( hirji_shifted<> & original_alg, double percent_tol = 1e-9 ) {

    // try a number of different n, smaller than the original n.
    BOOST_FOREACH( size_t factor, boost::counting_range( 1, 7 ) ) {

        // the smaller n
        const size_t n1 = original_alg.n * factor / 7;

        // is it worth testing?
        if( n1 > 0 && n1 != original_alg.n ) {

            //
            // Get the log p.m.f.
            //
            hirji_shifted< double >::slice_t log_pmf = original_alg.get_log_pmf( n1 );
            const size_t Q = boost::size( log_pmf );
            BOOST_CHECK_GT( Q, boost::math::round( n1 * ( - original_alg.min_log_b ) / original_alg.delta ) );
            cout
                << "Original distribution had n = " << original_alg.n << ", with Q = " << original_alg.get_Q() << ", "
                << "trying to recreate distribution for n1 = " << n1 << ", with Q = " << Q << "\n"
                ;

            //
            // Calculate the log c.m.f. for the original n
            //
            std::vector< double > original_cmf;
            {
                const double original_log_cumulative_total = make_cumulative( log_pmf, original_cmf );
                BOOST_CHECK_EQUAL( original_cmf.size(), std::vector< double >::size_type( Q ) );
                BOOST_CHECK_EQUAL( original_cmf.size(), std::vector< double >::size_type( boost::size( log_pmf ) ) );
                BOOST_CHECK_CLOSE( std::exp( original_log_cumulative_total ), 1., percent_tol );
            }

            //
            // Get the log p.m.f. using n1 < n
            //
            hirji_shifted<> smaller_n_alg( original_alg.b, n1, Q );

            //
            // Calculate the log c.m.f. for the new smaller n distribution
            //
            std::vector< double > smaller_cmf;
            {
                BOOST_CHECK_EQUAL( boost::size( smaller_n_alg.get_log_pmf( n1 ) ), int( Q ) );
                const double smaller_log_cumulative_total = make_cumulative( smaller_n_alg.get_log_pmf( n1 ), smaller_cmf );
                BOOST_CHECK_EQUAL( smaller_cmf.size(), Q );
                BOOST_CHECK_CLOSE( std::exp( smaller_log_cumulative_total ), 1., percent_tol );
            }

            //
            // Check the values are close
            //
            const double tol_lower_bound = boost::math::log1p( - .01 * percent_tol );
            const double tol_upper_bound = boost::math::log1p( + .01 * percent_tol );
            // cout << tol_lower_bound << " " << tol_upper_bound << "\n";
            for( int j = 1; int( Q - 1 ) >= j; ++j ) {
                // because of rounding errors we seem to be only able to get the c.m.f. to agree
                // up to +/- 2 lattice positions either side.
                // This effect appears to affect the unimportant end of the c.m.f. rather than
                // the end with low p-values.
                const int lower_j = std::max( 0, j - 2 );
                const int upper_j = std::min( int( Q - 1 ), j + 2 );
                //const double original_p = std::exp( original_cmf[ j ] );
                //const double smaller_n_p = std::exp( smaller_cmf[ j ] );
                // cout << j << "; original=" << original_cmf[ j ] << "; smaller=" << smaller_cmf[ j ] << "\n";
                BOOST_CHECK_GT( original_cmf[ j ], smaller_cmf[ lower_j ] + tol_lower_bound );
                BOOST_CHECK_LT( original_cmf[ j ], smaller_cmf[ upper_j ] + tol_upper_bound );
                BOOST_CHECK_LT( original_cmf[ lower_j ] + tol_lower_bound, smaller_cmf[ j ] );
                BOOST_CHECK_GT( original_cmf[ upper_j ] + tol_upper_bound, smaller_cmf[ j ] );
            }
        }
    }
}



BOOST_AUTO_TEST_CASE( log_add_test )
{
    BOOST_CHECK_CLOSE( log_add( std::log( 3. ), std::log( 9. ) ), std::log( 12. ), 1e-12 );
    BOOST_CHECK_CLOSE( log_add( std::log( 9. ), std::log( 3. ) ), std::log( 12. ), 1e-12 );
}


//BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( hirji_shifted_1, 2 )

BOOST_AUTO_TEST_CASE( hirji_shifted_1 )
{
    const double b[] = { .2, .2, .3, .3 };
    const std::vector< unsigned > ns = boost::assign::list_of( 1 )( 2 )( 5 )( 10 )( 20 );

    BOOST_FOREACH( unsigned n, ns ) {

        boost::timer t;
        hirji_shifted<> alg( b, n );
        cout << "Ours took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << n << " with Q=" << alg.get_Q() << "\n";

        check_is_proper_dist( alg, 1e-9 );
        check_for_smaller_n( alg );
    }
}


BOOST_AUTO_TEST_CASE( hirji_shifted_check_smaller )
{
    const double b[] = { .2, .2, .3, .3 };
    const std::vector< unsigned > ns = boost::assign::list_of( 20 );

    BOOST_FOREACH( unsigned n, ns ) {

        const size_t Q = 10 * n;
        boost::timer t;
        hirji_shifted<> alg( b, n, Q );
        cout << "Ours took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << n << " with Q=" << alg.get_Q() << "\n";

        check_is_proper_dist( alg, 1e-9 );
        check_for_smaller_n( alg );
    }
}


// do some longer tests for release builds
#ifdef NDEBUG

/// Check case where b's are all different
BOOST_AUTO_TEST_CASE( hirji_shifted_diff_bs )
{
    const double b[] = { .1, .2, .3, .4 };
    const std::vector< unsigned > ns = boost::assign::list_of( 1000 );

    BOOST_FOREACH( unsigned n, ns ) {

        const size_t Q = - std::log( .1 ) * 1e3 + 1;
        boost::timer t;
        hirji_shifted<> alg( b, n, Q );
        cout << "Ours took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << n << " with Q=" << alg.get_Q() << "\n";

        check_is_proper_dist( alg, 1e-9 );
        check_for_smaller_n( alg );
    }
}


BOOST_AUTO_TEST_CASE( hirji_shifted_large )
{
    const double b[] = { .2, .2, .3, .3 };
    //const double b[] = { .3, .3, .2, .2 };
    const unsigned n = 400;

    boost::timer t;
    hirji_shifted<> alg( b, n );
    cout << "Ours took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << n << " with Q=" << alg.get_Q() << "\n";
    check_is_proper_dist( alg, 1e-9 );
    check_for_smaller_n( alg );
}


BOOST_AUTO_TEST_CASE( hirji_shifted_extra_large )
{
    const double b[] = { .2, .2, .3, .3 };
    //const double b[] = { .3, .3, .2, .2 };
    const std::vector< unsigned > ns = boost::assign::list_of
        ( 1000 )
        // ( 2000 )
        // ( 3000 )
        // ( 5000 )
        // ( 7000 )
        // ( 10000 ) // memory requirements > 3Gb
        ;

    BOOST_FOREACH( unsigned n, ns ) {
        boost::timer t;
        hirji_shifted<> alg( b, n );
        cout << "Ours took " << t.elapsed() << " seconds to calculate p.m.f. for n=" << n << " with Q=" << alg.get_Q() << "\n";
        check_is_proper_dist( alg, 1.e-8 );
        check_for_smaller_n( alg, 1.e-8 );
    }
}


#endif //NDEBUG
