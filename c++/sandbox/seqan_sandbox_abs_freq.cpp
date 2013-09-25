/**
 * Copyright John Reid 2011
 */

#include <seqan/index.h>
#include <seqan/find_motif.h>
#include <iostream>
#include <stdexcept>

using namespace seqan;

typedef seqan::FrequencyDistribution< seqan::Dna5 >            freq_dist_t;        ///< The type that holds frequency distributions.
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
//			std::string tmp_meta;
//			assign( tmp_meta, meta );
//			ids.push_back( tmp_meta );
			appendValue( sequences, str );
		}
		std::cout << "Read " << length( sequences ) << " sequences with a total of " << num_bases << " bases\n";
		f.close();
	}

	index_t index( sequences );
	freq_dist_t seqan_dist;
	absFreqOfLettersInSetOfSeqs(
		seqan_dist,
		begin( getFibre( index, seqan::EsaText() ) ),
		end( getFibre( index, seqan::EsaText() ) )
	);
//	raw_text_t raw_text = getFibre( index, ESA_RawText() );
//	raw_text_it it = begin( raw_text );
//	++it;
	//--it; // compilation problem with operator--

	return 0;
}
