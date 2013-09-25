/**
 * Copyright John Reid 2010
 */

#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <vector>

template< typename Vec >
void
test_add_speed( Vec & vec, const size_t repetitions ) {
	typedef typename Vec::value_type value_t;

	// initialise vectors with values.
	for( size_t i = 0; vec.size() != i; ++i ) {
		vec[ i ] = double( i ) + 1./.3;
	}

	boost::timer timer;
	double sum = 0;
	for( size_t rep = 0; repetitions != rep; ++rep ) {
		BOOST_FOREACH( value_t v, vec ) {
			sum += v;
		}
	}
	const double elapsed = timer.elapsed();
	::std::cout << "Adding "<<repetitions*vec.size()<<" took "<<elapsed<<" seconds\n";
}


int
main( int argc, char * argv[] ) {

	typedef std::vector< int >        int_vec;
	typedef std::vector< double >     double_vec;

	const size_t vec_size = 5000;
	const size_t repetitions = 1e6;

	int_vec ints( vec_size );
	double_vec doubles( vec_size );

	::std::cout << "Adding integers...\n";
	test_add_speed( ints, repetitions );
	::std::cout << "Adding doubles...\n";
	test_add_speed( doubles, repetitions );

	return 0;
}
