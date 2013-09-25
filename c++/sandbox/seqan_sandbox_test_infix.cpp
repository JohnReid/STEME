/**
 * Copyright John Reid 2011
 *
 * \file
 * Code to test infix of Concatenator
 */

#include <seqan/sequence.h>
#include <stdexcept>

using namespace seqan;

typedef String< Dna5 >                        string_t;
typedef StringSet< string_t >                 string_set_t;
typedef ConcatenatorManyToOne< string_set_t > concatenator_t;
typedef Infix< concatenator_t >::Type         infix_t;

int
main( int argc, char * argv [] ) {
	string_set_t seqs;
	appendValue( seqs, string_t( "ACGT" ) );
	concatenator_t concat( seqs );
	infix_t in( concat, 0, 4 );
	if( 4 != length( in ) ) throw std::logic_error( "length() not working." );
	if( difference( begin( in ), end( in ) ) != int( length( in ) ) ) throw std::logic_error( "difference() not working." );
	return 0;
}
