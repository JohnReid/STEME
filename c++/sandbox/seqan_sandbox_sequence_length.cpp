/**
 * Copyright John Reid 2010
 */

#include <seqan/index.h>

int
main( int argc, char * argv[] )
{
	using namespace seqan;

	typedef String< Dna5 >                string_t;        ///< A string of Dna5.
	typedef seqan::StringSet< string_t >  string_set_t;    ///< StringSet type.

	sequenceLength( 0, string_set_t() );

	return 0;
}
