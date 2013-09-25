/**
 * Copyright John Reid 2011
 */

#include <seqan/index.h>
#include <boost/assert.hpp>
#include <iostream>

using namespace seqan;

typedef String< Dna5 >                                                string_t;              /**< A string of Dna5. */
typedef seqan::StringSet< string_t >                                  string_set_t;          /**< StringSet type. */
typedef Index< string_set_t >                                         index_t;               /**< The index over our strings. */
typedef seqan::Iterator< index_t, seqan::TopDown<> >::Type            top_down_it;           /**< A iterator over the index type. */




int
main( int argc, char * argv[] ) {

	string_set_t sequences;
	appendValue( sequences, "ACGACTACGAGC" );
	appendValue( sequences, "ACCCAAC" );
	appendValue( sequences, "" );
	index_t index( sequences );

	String< double > property;
	const double v = 1.3;
	resizeVertexMap( index, property );
	assignProperty( property, value( top_down_it( index ) ), v );
	const double prop_value = getProperty( property, value( top_down_it( index ) ) );
	std::cout << v << " == " << prop_value << "\n";
	BOOST_ASSERT( v == prop_value );

	return 0;
}
