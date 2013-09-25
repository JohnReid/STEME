/** Copyright John Reid 2011
 *
 * \file
 * \brief Test compilation of steme header.
 */



#include <steme/find_best_w_mers.h>
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;

void
empty_test()
{
}

test_suite *
init_unit_test_suite( int argc, char * argv [] )
{
	framework::master_test_suite().add( BOOST_TEST_CASE( &empty_test ) );
	return 0;
}


