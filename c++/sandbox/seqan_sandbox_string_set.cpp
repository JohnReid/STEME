/**
 * Copyright John Reid 2011
 */

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <fstream>

using namespace seqan;

typedef String< Dna5 >                                         string_t;              /**< A string of Dna5. */
typedef StringSet< string_t >                                  string_set_t;          /**< StringSet type. */


/**
 * Read a FASTA file into a string set.
 */
void
read_fasta( const char * filename, string_set_t & sequences )
{
	::std::ifstream f( filename, ::std::ios_base::in | ::std::ios_base::binary );
	if( ! f ) {
		std::cerr << "Could not open FASTA file: " << filename << "\n";
	} else {
		String< char > meta;
		string_t str;
		while( f && f.peek() != ::std::char_traits< char >::eof() ) {
			readMeta( f, meta, Fasta() );
			read( f, str, Fasta() );
			appendValue( sequences, str );
		}
	}
	f.close();
}


void show_ids( const string_set_t & sequences ) {
	int pos = 0;
	for( Iterator< const string_set_t >::Type i = begin( sequences ); end( sequences ) != i; ++i, ++pos ) {
		std::cout << positionToId( sequences, pos ) << "\n";
		//std::cout << getValue( i ) << "\n";
	}
}


int
main( int argc, char * argv[] ) {

	if( argc < 2 ) {
		std::cerr << "USAGE: " << argv[0] << " <fasta file>\n";
		return -1;
	} else {
		string_set_t sequences;
		read_fasta( argv[1], sequences );
		removeValueById( sequences, 2 );
		appendValue( sequences, string_t( "ACGTN" ) );
		show_ids( sequences );
		std::cout << getValueById( sequences, 3 ) << "\n";
		return 0;
	}
}
