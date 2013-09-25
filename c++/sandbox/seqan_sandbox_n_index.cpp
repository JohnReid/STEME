/**
 * Copyright John Reid 2009-2011
 */

#include <seqan/index.h>

using namespace seqan;

typedef String< Dna5 >                             string_t;              /**< A string of Dna5. */
typedef StringSet< string_t >                      string_set_t;          /**< StringSet type. */
typedef Index< string_set_t >                      index_t;               /**< The index over our strings. */
typedef Finder< index_t >                          finder_t;              /**< Finder type. */
typedef Iterator< index_t, TopDown<> >::Type       top_down_it;           /**< A iterator over the index type. */
typedef Position< finder_t >::Type                 position_t;            ///< A position in the strings.

int
main( int argc, char * argv[] ) {

	string_set_t sequences;
	appendValue( sequences, "AACNACGT" );
	appendValue( sequences, "TGAACTG" );
	appendValue( sequences, "ANACAAC" );
	index_t index( sequences );
	finder_t finder( index );
	while( find( finder, "NAC" ) ) {
		position_t pos = position( finder );
	    ::std::cout << pos << " : " << getSeqNo( pos ) << ", " << getSeqOffset( pos ) << "\n";
	    ::std::cout << posGlobalize( pos, stringSetLimits( sequences ) ) << "\n";
	}

	top_down_it it( index );
	goDown( it );
	goDown( it );
	::std::cout << prefix( representative( it ), 3 ) << "\n";
	typedef Infix< Fibre< index_t, EsaSA >::Type const >::Type occurrences_t;
	typedef Value< occurrences_t >::Type occurrence_t;
	const occurrences_t occs = seqan::getOccurrences( it );
	for( unsigned i = 0; seqan::length( occs ) != i; ++i ) {
		const occurrence_t & occ = occs[ i ];
		::std::cout << getSeqNo( occ ) << ", " << getSeqOffset( occ ) << "\n";
	}

	return 0;
}
