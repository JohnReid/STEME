/**
 * Copyright John Reid 2009
 */

#include <seqan/index.h>

using namespace seqan;

typedef String< Dna5 >                                        string_t;              /**< A string of Dna5. */
typedef seqan::StringSet< string_t >                          string_set_t;          /**< StringSet type. */
typedef Index< string_set_t >                                 index_t;               /**< The index over our strings. */
typedef Fibre< index_t, EsaRawText >::Type                    raw_text_t;            /**< The type of the raw text (concatentation of the text) in the index. */
typedef Iterator< raw_text_t >::Type                          raw_text_it;           /**< The type of the raw text iterator. */

int
main( int argc, char * argv[] ) {

	string_set_t sequences;
	index_t index( sequences );
	raw_text_t raw_text = getFibre( index, EsaRawText() );
	raw_text_it it = begin( raw_text );
	++it;
	//--it; // compilation problem with operator--

	return 0;
}
