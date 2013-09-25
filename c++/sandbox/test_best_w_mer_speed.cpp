/**
 * Copyright John Reid 2010
 */

#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <iterator>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
//#include "sti/sset.h"

typedef boost::mt19937                                rng_t;     // RNG type.
typedef boost::uniform_01<>                           dist_t;    // Distribution type.
typedef boost::variate_generator< rng_t &, dist_t >   gen_t;     // Generator type.

rng_t rng( 42u );
dist_t dist;
gen_t gen( rng, dist );

const size_t num_to_save = 200;
const size_t num_to_generate = 1000000;
const size_t num_repetitions = 100;

typedef std::vector< double > data_vec;
data_vec data;

void generate_data() {
	data.resize( num_to_generate );
	for( size_t i = 0; num_to_generate != i; ++i ) {
		data[ i ] = gen();
	}
}

template< typename set = std::set< double > >
struct best_set {
	set best;

	best_set()
	{ }

	double worst() {
		return best.empty()
			? -std::numeric_limits< double >::max()
			: *best.begin();
	}
	size_t size() { return best.size(); }

	void insert( double x ) {
		best.insert( x );
		if( size() > num_to_save ) {
			best.erase( best.begin() );
		}
	}
};


template< typename T >
struct greater_than {
	bool operator()( T x1, T x2 ) { return x1 > x2; }
};


template< typename vec = std::vector< double > >
struct best_vector {
	vec best;

	best_vector()
	{
		best.reserve( num_to_save + 1 );
	}

	double worst() {
		return best.empty()
			? -std::numeric_limits< double >::max()
			: *best.rbegin();
	}
	size_t size() { return best.size(); }

	void insert( double x ) {
		best.insert(
		    std::upper_bound( best.begin(), best.end(), x, greater_than< double >() ),
		    x
		);

		if( size() > num_to_save ) {
			best.resize( num_to_save );
		}
	}
};



template< typename Container >
Container test_insertion() {

	Container s;
	double sum = 0.;
	for( size_t i = 0; num_to_generate != i; ++i ) {

		const double x = data[ i ];

		// work out where we would insert this data
		if( s.size() < num_to_save || x > s.worst() ) {
			s.insert( x );
			BOOST_FOREACH( double y, s.best ) {
				sum += y;
			}
		}
	}

	return s;
}

template< typename Container >
Container test_repeat_insertion() {
	Container result;
	for( size_t i = 0; num_repetitions != i; ++i ) {
		result = test_insertion< Container >();
	}
	return result;
}






int
main( int argc, char * argv[] ) {

	generate_data();

	boost::timer timer;
	double elapsed;

	timer.restart();
	best_vector<> vec = test_repeat_insertion< best_vector<> >();
	elapsed = timer.elapsed();
	::std::cout << "Adding "<<num_to_generate<<" to a vector "<<num_repetitions<<" times took "<<elapsed<<" seconds\n";
	//std::copy( vec.best.begin(), vec.best.end(), std::ostream_iterator< double >( std::cout, "\n" ) );

	timer.restart();
	best_set<> set = test_repeat_insertion< best_set<> >();
	elapsed = timer.elapsed();
	::std::cout << "Adding "<<num_to_generate<<" to a set    "<<num_repetitions<<" times took "<<elapsed<<" seconds\n";
	//std::copy( set.best.begin(), set.best.end(), std::ostream_iterator< double >( std::cout, "\n" ) );

//	timer.restart();
//	best_set< sti::set< double > > sti_set = test_repeat_insertion< best_set< sti::set< double > > >();
//	elapsed = timer.elapsed();
//	::std::cout << "Adding "<<num_to_generate<<" to a set    "<<num_repetitions<<" times took "<<elapsed<<" seconds\n";
	//std::copy( set.best.begin(), set.best.end(), std::ostream_iterator< double >( std::cout, "\n" ) );

	return 0;
}
