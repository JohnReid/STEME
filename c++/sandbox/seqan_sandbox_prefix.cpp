/**
 * Copyright John Reid 2009-2011
 */

#include <seqan/index.h>

using namespace seqan;

typedef String< Dna5 >                                                string_t;              /**< A string of Dna5. */
typedef seqan::StringSet< string_t >                                  string_set_t;          /**< StringSet type. */
typedef Index< string_set_t >                                         index_t;               /**< The index over our strings. */
typedef seqan::Iterator< index_t, seqan::TopDown<> >::Type            top_down_it;           /**< A iterator over the index type. */
typedef Fibre< index_t, EsaText >::Type                               text_t;                /**< The type of the text in the index. */
typedef Infix< text_t const >::Type                                   representative_t;      /**< Representative type. */
typedef Prefix< representative_t >::Type                              prefix_t;              /**< Type of a prefix of a representative. */




int
main( int argc, char * argv[] ) {

	string_set_t sequences;
	index_t index( sequences );
	top_down_it it( index );
	representative_t repr = representative( it );
	prefix_t p = prefix( repr, 0 );

	return 0;
}
