/**
 * Copyright John Reid 2011
 */

#include <seqan/index.h>

#include <boost/assert.hpp>
#include <boost/tuple/tuple.hpp>

#include <iostream>
#include <limits>

using namespace seqan;

typedef String< Dna5 >                                  string_t;              /**< A string of Dna5. */
typedef StringSet< string_t >                           string_set_t;          /**< StringSet type. */
typedef Index< string_set_t >                           index_t;               /**< The index over our strings. */
typedef Iterator< index_t, TopDown<> >::Type            top_down_it;           /**< A iterator over the index type. */
typedef boost::tuple< int, int >                        range_t;               ///< Range of ints.
typedef SAValue< index_t >::Type                        sa_value_t;            ///< Suffix array default alphabet type.
typedef Fibre< index_t, EsaSA >::Type                   fibre_sa_t;            ///< Fibre_SA type.
typedef Fibre< index_t, EsaText >::Type                 text_t;                ///< Fibre_SA type.
typedef Infix< fibre_sa_t const >::Type                 occurrences_t;         ///< Stores a sequence of positions.
typedef Infix< text_t const >::Type                     representative_t;      ///< Representative type.
typedef Prefix< representative_t >::Type                prefix_t;              ///< Type of a prefix of a representative.

void
update_range( range_t & range, int pos ) {
	if( pos < range.get< 0 >() ) {
		range.get< 0 >() = pos;
	}
	if( pos > range.get< 1 >() ) {
		range.get< 1 >() = pos;
	}
}


void
visit_node( top_down_it it, range_t & range ) {

	std::cout << representative( it );
	for( size_t i = 0; 15 - repLength( it ) != i; ++i ) std::cout << " ";
	occurrences_t occs = seqan::getOccurrences( it );
	typedef Size< occurrences_t >::Type occ_size_t;
	const occ_size_t num_occurrences = seqan::length( occs );
	for( occ_size_t i = 0; num_occurrences != i; ++i ) {
		const int global_pos = seqan::posGlobalize( occs[ i ], stringSetLimits( getFibre( container( it ), EsaText() ) ) );
		update_range( range, global_pos );
		std::cout << ", " << global_pos;
		std::cout << ":" << getValueI1( occs[ i ] ) << ":" << getValueI2( occs[ i ] );
	}
	std::cout << "\n";

	if( goDown( it ) ) {
		visit_node( it, range );
		while( goRight( it ) ) {
			visit_node( it, range );
		}
	}
}



int
main( int argc, char * argv[] ) {

	string_set_t sequences;
	appendValue( sequences, "ACGACTACGAGC" );
	appendValue( sequences, "ACCCAAC" );
	appendValue( sequences, "" );
	index_t index( sequences );

	range_t range(
		std::numeric_limits< int >::max(),
		-std::numeric_limits< int >::max()
	);
	visit_node( top_down_it( index ), range );
	std::cout << range.get< 0 >() << " - " << range.get< 1 >() << "\n";

	return 0;
}
