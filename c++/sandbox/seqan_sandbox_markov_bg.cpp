/**
 * Copyright John Reid 2011
 */

#include <seqan/index.h>
#include <seqan/statistics.h>
#include <iostream>
#include <stdexcept>

using namespace seqan;

typedef String< Dna5 >                                        string_t;              /**< A string of Dna5. */
typedef seqan::StringSet< string_t >                          string_set_t;          /**< StringSet type. */
typedef Index< string_set_t >                                 index_t;               /**< The index over our strings. */
typedef Fibre< index_t, EsaRawText >::Type                    raw_text_t;            /**< The type of the raw text (concatentation of the text) in the index. */
typedef Iterator< raw_text_t >::Type                          raw_text_it;           /**< The type of the raw text iterator. */

#define FASTA "/home/john/Dev/MyProjects/STEME/python/test/fasta/T00759-small.fa"

int
main( int argc, char * argv[] ) {

	string_set_t sequences;

	::std::ifstream f( FASTA, ::std::ios_base::in | ::std::ios_base::binary );
	size_t num_bases = 0;
	if( ! f ) {
		throw std::logic_error( "Could not open: \"" FASTA "\"" );
	} else {
		String< char > meta;
		string_t str;
		while( f ) {
			readMeta( f, meta, Fasta() );
			read( f, str, Fasta() );
			if( length( str ) == 0 && length( meta ) == 0 ) continue;
			num_bases += length( str );
			appendValue( sequences, str );
		}
		std::cout << "Read " << length( sequences ) << " sequences with a total of " << num_bases << " bases\n";
		f.close();
	}

	index_t index( sequences );
	const size_t order = 0;
	MarkovModel< Dna5 > mm( order );
	mm.build( sequences );
	string_t s = "N";
	std::cout << "p(" << s << ")=" << mm.emittedProbability( s ) << "\n";
	s = "A";
	std::cout << "p(" << s << ")=" << mm.emittedProbability( s ) << "\n";
	s = "C";
	std::cout << "p(" << s << ")=" << mm.emittedProbability( s ) << "\n";
	s = "G";
	std::cout << "p(" << s << ")=" << mm.emittedProbability( s ) << "\n";
	s = "T";
	std::cout << "p(" << s << ")=" << mm.emittedProbability( s ) << "\n";

	return 0;
}
