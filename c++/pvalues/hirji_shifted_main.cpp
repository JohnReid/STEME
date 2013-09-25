/**
 * Copyright John Reid 2011
 *
 * \file Code that calculates the p.m.f. and c.m.f. using the Hirji shifted algorithm.
 */

#include <pvalues/hirji_shifted.h>

#include <boost/timer.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <vector>
#include <iostream>
#include <numeric>
#include <algorithm>

using namespace boost;
using namespace std;


int
main( int argc, const char * argv[] ) {

    if( argc < 5 ) {
        cerr
            << "USAGE: " << argv[ 0 ]
            << " <N> <Q> <b_1> ... <b_A>\n"
            ;
        return -1;
    }

    const size_t N = lexical_cast< size_t >( argv[ 1 ] );
    size_t Q = lexical_cast< size_t >( argv[ 2 ] );
    double( * arg_to_double )( const std::string & ) = lexical_cast< double >;
    std::vector< double > bg(
        make_transform_iterator( argv + 3   , arg_to_double ),
        make_transform_iterator( argv + argc, arg_to_double )
    );
    for_each( bg, boost::lambda::_1 = boost::lambda::_1 / std::accumulate( bg.begin(), bg.end(), 0. ) );

    if( ! Q ) {
        Q = 2 * N;
    }
    cout << "N = " << N << "\n";
    cout << "Q = " << Q << "\n";
    cout << "background = ";
    copy( bg.begin(), bg.end(), ostream_iterator< double >( cout, " " ) );
    cout << "\n";

    timer t;
    pvalues::hirji_shifted<> alg( bg, N, Q );
    const double elapsed = t.elapsed();
    cout << "delta = " << alg.delta << "\n";
    const pvalues::hirji_shifted<>::slice_t log_pmf = alg.get_log_pmf( N );
    for( size_t j = 0; Q != j; ++j ) {
        cout << "log p(round(LLR/delta)=" << j << ") = " << log_pmf[ j ] << "\n";
    }
    cout << "Took " << elapsed << " seconds to calculate probability mass function.\n";
}



