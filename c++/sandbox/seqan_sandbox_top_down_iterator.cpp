/**
 * Copyright John Reid 2009, 2010, 2011
 */

#include <seqan/index.h>
#include <seqan/sequence.h>

using namespace seqan;

typedef String< Dna5 >                                                string_t;              /**< A string of Dna5. */
typedef seqan::StringSet< string_t >                                  string_set_t;          /**< StringSet type. */
typedef Index< string_set_t >                                         index_t;               /**< The index over our strings. */
typedef seqan::Iterator< index_t, seqan::TopDown<> >::Type            top_down_it;           /**< A iterator over the index type. */


/**
 * Visit a node in the suffix tree (array).
 */
void
visit( top_down_it it, int last_depth ) {
	std::cout <<  std::setw( 6 ) << countOccurrences( it ) << " occurrences of \"" << representative( it ) << "\"\n";
	const int new_depth = repLength( it );
	assert( new_depth > last_depth );
	if( goDown( it ) ) {
		std::cout << "Gone down\n";
		visit( it, new_depth );
		while( goRight( it ) ) {
			std::cout << "Gone right\n";
			visit( it, new_depth );
		}
	}
}


/**
 * Read a FASTA file into a string set.
 */
void
read_fasta( const char * filename, string_set_t & sequences )
{
	::std::ifstream f( filename );
	if( ! f ) {
		std::cerr << "Could not open FASTA file: " << filename << "\n";
	} else {
		String< char > meta;
		string_t str;
		while( f ) {
			readMeta( f, meta, Fasta() );
			read( f, str, Fasta() );
			appendValue( sequences, str );
		}
	}
	f.close();
}


int
main( int argc, char * argv[] ) {

	if( argc < 2 ) {
		std::cerr << "USAGE: " << argv[0] << " <fasta file>\n";
		return -1;
	} else {
		string_set_t sequences;
		read_fasta( argv[1], sequences );
		index_t index( sequences );
		visit( top_down_it( index ), -1 );
		return 0;
	}
}
