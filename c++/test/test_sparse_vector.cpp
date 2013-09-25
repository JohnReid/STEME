/** Copyright John Reid 2011, 2012, 2013
 *
 * \file
 * \brief Test sparse vector iteration.
 */



#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <numeric>
#include <functional>
#include <cmath>
#include <algorithm>

using namespace boost::unit_test;

template< typename SparseVector >
void
write_iterator( SparseVector & vec, typename SparseVector::iterator i ) {
	if( vec.end() == i ) {
		std::cout << "<END>";
	} else {
		std::cout << i.index() << " : " << *i;
	}
}

template< typename SparseVector >
void
test_sparse_vector()
{
	typedef typename SparseVector::iterator iterator;
	typedef typename SparseVector::value_type value_type;

	const size_t size = 10;

	SparseVector sparse_vec( size );
	std::vector< value_type > dense_vec( size );

	sparse_vec[ 3 ] = dense_vec[ 3 ] = 1;
	sparse_vec[ 8 ] = dense_vec[ 8 ] = 2;

	// check basic iteration
	for( iterator i = sparse_vec.begin(); sparse_vec.end() != i; ++i ) {
        write_iterator< SparseVector >( sparse_vec, i );
		std::cout << "\n";
	}
	BOOST_CHECK_EQUAL( 1, *sparse_vec.begin() );
	BOOST_CHECK_EQUAL( 2, *(++sparse_vec.begin()) );

	// check find...
	for( size_t i = 0; size != i; ++i ) {
		std::cout << "find(" << i << ") - ";
		iterator it = sparse_vec.find( i );
		write_iterator< SparseVector >( sparse_vec, it );
		std::cout << "\n";
	}
	BOOST_CHECK_EQUAL( 1, *sparse_vec.find( 0 ) );
	BOOST_CHECK_EQUAL( 1, *sparse_vec.find( 1 ) );
	BOOST_CHECK_EQUAL( 1, *sparse_vec.find( 2 ) );
	BOOST_CHECK_EQUAL( 1, *sparse_vec.find( 3 ) );
	BOOST_CHECK_EQUAL( 2, *sparse_vec.find( 4 ) );
	BOOST_CHECK_EQUAL( 2, *sparse_vec.find( 5 ) );
	BOOST_CHECK_EQUAL( 2, *sparse_vec.find( 6 ) );
	BOOST_CHECK_EQUAL( 2, *sparse_vec.find( 7 ) );
	BOOST_CHECK_EQUAL( 2, *sparse_vec.find( 8 ) );
	BOOST_CHECK( sparse_vec.end() == sparse_vec.find( 9 ) );

	// check accumulate...
	BOOST_CHECK_EQUAL( 0, std::accumulate( sparse_vec.find( 3 ), sparse_vec.find( 3 ), 0 ) );
	BOOST_CHECK_EQUAL( 1, std::accumulate( sparse_vec.find( 3 ), sparse_vec.find( 8 ), 0 ) );
	BOOST_CHECK_EQUAL( 3, std::accumulate( sparse_vec.find( 3 ), sparse_vec.find( 9 ), 0 ) );
	BOOST_CHECK_EQUAL( 3, std::accumulate( sparse_vec.find( 0 ), sparse_vec.find( 9 ), 0 ) );

	// check max...
	const value_type &( * maxfn )( const value_type &, const value_type & ) = std::max< value_type >;
	BOOST_CHECK_EQUAL( 0, std::accumulate( sparse_vec.find( 3 ), sparse_vec.find( 3 ), 0, maxfn ) );
	BOOST_CHECK_EQUAL( 1, std::accumulate( sparse_vec.find( 3 ), sparse_vec.find( 8 ), 0, maxfn ) );
	BOOST_CHECK_EQUAL( 2, std::accumulate( sparse_vec.find( 3 ), sparse_vec.find( 9 ), 0, maxfn ) );
	BOOST_CHECK_EQUAL( 2, std::accumulate( sparse_vec.find( 0 ), sparse_vec.find( 9 ), 0, maxfn ) );
}

test_suite *
init_unit_test_suite( int argc, char * argv [] )
{
	framework::master_test_suite().add( BOOST_TEST_CASE( &test_sparse_vector< boost::numeric::ublas::mapped_vector< int > > ) );
	framework::master_test_suite().add( BOOST_TEST_CASE( &test_sparse_vector< boost::numeric::ublas::compressed_vector< int > > ) );
	// cannot compile this one with the sparse_vec.find() call.
	//framework::master_test_suite().add( BOOST_TEST_CASE( &test_sparse_vector< boost::numeric::ublas::coordinate_vector< int > > ) );
	return 0;
}

