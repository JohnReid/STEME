/**
 * Copyright John Reid 2010
 */

#include <seqan/file.h>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

int
main( int argc, char * argv[] ) {

	using namespace seqan;

	const char * fasta_filename = "/home/john/Dev/MyProjects/Bio/MotifSearch/python/test/fasta/find-best-w-mers-test.fa";
	::std::cout << "Opening " << fasta_filename << "\n";
	::std::ifstream fstrm( fasta_filename, ::std::ios_base::in | ::std::ios_base::binary );
	if( ! fstrm ) {
		throw ::std::logic_error( "Could not open fasta file" );
	}
	String< Dna > seq;
	String< char > meta;
	size_t num_bases = 0;
	size_t num_seqs = 0;
	while( ! fstrm.eof() )
	{
		readMeta( fstrm, meta, Fasta() );
		read( fstrm, seq, Fasta() );
		if( length( seq ) == 0 && length( meta ) == 0 ) continue;
		num_bases += length( seq );
		::std::cout << "Read " << length( seq ) << " bases into sequence \"" << meta << "\"\n";
		++num_seqs;
	}
	::std::cout << "Read " << num_seqs << " sequences with a total of " << num_bases << " bases\n";


	return 0;
}
