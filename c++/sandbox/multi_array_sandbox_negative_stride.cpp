/**
 * Copyright John Reid 2009
 */

//#define BOOST_DISABLE_ASSERTS

#include <boost/multi_array.hpp>
#include <iostream>

int
main( int argc, char * argv[] )
{
	using namespace boost;

	const size_t size = 4;

	typedef multi_array< double, 1 >            array_type;
	array_type my_array( extents[ size ] );
	my_array[0] = 1.;
	my_array[1] = 2.;
	my_array[2] = 3.;
	my_array[3] = 4.;

	typedef multi_array_types::index_range      range;
	array_type::array_view< 1 >::type my_view = my_array[ indices[ range( 3, -1, -1 ) ] ];
	//array_type::array_view< 1 >::type my_view = my_array[ indices[ range().start( 3 ).stride( -1 ) ] ];

	for( size_t i = 0; my_view.size() != i; ++i ) {
		std::cout << i << " : " << my_view[ i ] << "\n";
	}

	return 0;
}
