/**
 * Copyright John Reid 2009
 */

#include <seqan/index.h>

using namespace seqan;

typedef String< Dna5 >                                 string_t;              ///< A string of Dna5.
typedef StringSet< string_t >                          string_set_t;          ///< StringSet type.
typedef Reference< string_set_t >::Type                string_ref_t;          ///< Reference to a string in a string set.
typedef Infix< string_t >::Type                        ref_infix_t;           ///< Infix of a string

int
main( int argc, char * argv[] ) {

	string_set_t sequences;
	string_ref_t string = getValueById( sequences, positionToId( sequences, 0 ) );
	ref_infix_t my_infix = infix( string, 3, 6 );

	return 0;
}
