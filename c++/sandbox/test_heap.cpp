/** Copyright John Reid 2009, 2010, 2011
 *
 * \file
 * \brief Tests heap storage for find best W-mers.
 *
 */

#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/timer.hpp>
#include <boost/test/utils/wrap_stringstream.hpp>

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <set>
#include <limits>

#define MAKE_STRING( x ) ( boost::wrap_stringstream().ref() << x ).str()

typedef std::vector< double > Z_vec;

void
generate_random_Z( Z_vec & Zs, size_t num, uint32_t seed = 1 ) {
	boost::random::mt19937 rng( seed );
	boost::random::uniform_01<> uniform;
	Zs.clear();
	Zs.reserve( num );
	for( size_t i = 0; num != i; ++i ) {
		Zs.push_back( uniform( rng ) );
	}
}


/**
 * Stores the best N elements added to it. Uses a heap based on vector storage. Better implementation below.
 */
template< typename T >
struct best_elements_heap_inefficient {
	typedef std::vector< T > storage_t;
	typedef typename storage_t::iterator iterator;

	size_t N;
	storage_t storage;

	template< typename InputIterator >
	best_elements_heap_inefficient( size_t N, InputIterator begin, InputIterator end )
	: N( N )
	{
		storage.reserve( N + 1 );
		storage.push_back( -std::numeric_limits< T >::max() );
		storage.insert( storage.end(), begin, end );
		BOOST_ASSERT( storage.size() <= N + 1 );
		heapify();
	}

	iterator
	begin() {
		return storage.begin() + 1;
	}

	iterator
	end() {
		return storage.end();
	}

	size_t
	size() {
		return end() - begin();
	}

	void
	heapify() {
		std::make_heap( storage.begin(), storage.end(), std::greater< T >() );
	}

	void
	add( T t ) {
		if( storage[0] < t ) {
			storage[0] = t;
			heapify();
		}
	}

	void
	sort() {
		std::sort( begin(), end() );
	}
};


/**
 * Stores the best N elements added to it. Uses a heap based on vector storage.
 */
template< typename T, typename _Compare = std::less< T > >
struct best_elements_heap {
	typedef std::vector< T > storage_t;
	typedef typename storage_t::reverse_iterator iterator;

	size_t N;
	storage_t storage;
	std::binary_negate< _Compare > compare;

	/// Constructor
	best_elements_heap( size_t N, _Compare compare = _Compare() )
	: N( N )
	, compare( std::not2( compare ) )
	{
		storage.reserve( N );
	}

	/// More efficient constructor that will take up to N values
	template< typename InputIterator >
	best_elements_heap( size_t N, InputIterator begin, InputIterator end, _Compare compare = _Compare() )
	: N( N )
	, compare( std::not2( compare ) )
	{
		storage.reserve( N );
		storage.insert( storage.end(), begin, end );
		if( storage.size() > N ) {
			throw std::logic_error( "Can only insert up to N values in constructor." );
		}
		std::make_heap( storage.begin(), storage.end(), compare );
	}

	inline
	iterator
	begin() {
		return storage.rbegin();
	}

	inline
	iterator
	end() {
		return storage.rend();
	}

	inline
	size_t
	size() {
		return storage.size();
	}

	inline
	T
	worst() const {
		return storage[ 0 ];
	}

	inline
	void
	add( T t ) {
		if( size() < N ) {
			storage.push_back( t );
			std::push_heap( storage.begin(), storage.end(), compare );
		} else {
			BOOST_ASSERT( N == size() );
			if( compare( t, worst() ) ) {
				std::pop_heap( storage.begin(), storage.end(), compare );
				*storage.rbegin() = t;
				std::push_heap( storage.begin(), storage.end(), compare );
			}
		}
	}

	inline
	void
	sort() {
		std::sort_heap( storage.begin(), storage.end(), compare );
	}
};


/**
 * Stores the best N elements added to it. Uses a set.
 */
template< typename T, typename _Compare = std::less< T > >
struct best_elements_set {
	typedef std::multiset< T, _Compare > storage_t;
	typedef typename storage_t::iterator iterator;

	size_t N;
	storage_t storage;

	template< typename InputIterator >
	best_elements_set( size_t N, InputIterator begin, InputIterator end )
	: N( N )
	, storage( begin, end )
	{
		BOOST_ASSERT( storage.size() <= N );
	}

	size_t
	size() {
		return storage.size();
	}

	iterator
	begin() {
		return storage.begin();
	}

	iterator
	end() {
		return storage.end();
	}

	void
	add( T t ) {
		BOOST_ASSERT( storage.size() == N );
		iterator worst = begin();
		if( *worst < t ) {
			storage.erase( worst );
			storage.insert( t );
#ifndef NDEBUG
			check_set();
#endif //NDEBUG
		}
		BOOST_ASSERT( storage.size() == N );
	}

	void
	check_set() {
		iterator i = begin();
		const double worst = *i;
		for( ++i; end() != i; ++i ) {
			BOOST_ASSERT( *i >= worst );
		}
	}

	void
	sort() {
	}
};


template< typename BE >
void
test_best_elements( const char * name, const Z_vec & Zs, const Z_vec & sorted_Zs, size_t max_size, bool show = false )
{
	// time the method
	boost::timer timer;

	// create the method object and populate with Zs
	Z_vec::const_iterator j = sorted_Zs.begin() + max_size;
	BE best( max_size, sorted_Zs.begin(), j );
	for( ; sorted_Zs.end() != j; ++j ) {
		best.add( *j );
	}
	std::cout << name << " took " << timer.elapsed() << " seconds\n";

	// sort best elements
	best.sort();

	// show if requested
	if( show ) {
		std::copy( best.begin(), best.end(), std::ostream_iterator< double >( std::cout, "," ) );
		std::cout << "\n";
	}

	// check right number of elements
	if( best.size() != max_size ) {
		throw std::logic_error( "Got wrong number of best elements" );
	}

	// check best elements match sorted ones.
	typename BE::iterator i = best.begin();
	BOOST_ASSERT( sorted_Zs.end() - sorted_Zs.begin() >= int( max_size ) );
	j = sorted_Zs.end() - max_size;
	while( best.end() != i ) {
		if( *i != *j ) {
			throw std::logic_error( MAKE_STRING( name << ": did not calculate best elements correctly: " << *i << " != " << *j ) );
		}
		++i;
		++j;
	}
	BOOST_ASSERT( sorted_Zs.end() == j );
}


int
main( int argc, char * argv[] ) {

#ifndef NDEBUG
	const size_t max_size = 20;
	const size_t num_Zs = 100;
	const bool show = true;
#else //NDEBUG
	const size_t max_size = 20;
	const size_t num_Zs = 1e7;
	const bool show = false;
#endif //NDEBUG

	std::cout << "Generating " << num_Zs << " random Zs\n";
	Z_vec Zs;
	generate_random_Z( Zs, num_Zs );

	std::cout << "Sorting " << num_Zs << " Zs\n";
	Z_vec sorted_Zs;
	{
		boost::timer timer;
		sorted_Zs.reserve( Zs.size() );
		std::copy( Zs.begin(), Zs.end(), std::back_insert_iterator< Z_vec >( sorted_Zs ) );
		std::sort( sorted_Zs.begin(), sorted_Zs.end() );
		std::cout << "sorting took " << timer.elapsed() << " seconds\n";
		if( show ) {
			std::copy( sorted_Zs.end() - max_size, sorted_Zs.end(), std::ostream_iterator< double >( std::cout, "," ) );
			std::cout << "\n";
		}
	}

	std::cout << "Looking for " << max_size << " best elements.\n";
	test_best_elements< best_elements_heap< double > >( "best_elements_heap", Zs, sorted_Zs, max_size, show );
	test_best_elements< best_elements_set< double > >( "best_elements_set", Zs, sorted_Zs, max_size, show );
	//test_best_elements< best_elements_heap_inefficient< double > >( "best_elements_heap_inefficient", Zs, sorted_Zs, max_size, show );

}
