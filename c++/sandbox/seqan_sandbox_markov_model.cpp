/**
 * Copyright John Reid 2009-2011
 */

#include <seqan/index.h>
#include <seqan/statistics.h>

using namespace seqan;

int main()
{
	const size_t order = 0;
	MarkovModel< Dna5 > mm( order );
	typedef String< Dna5 > string_t;
	StringSet< string_t > X;
	appendValue( X, string_t( "AATTAAT" ) );
	appendValue( X, string_t( "CGGGCGGGGCGCC" ) );
	appendValue( X, string_t( "CGATCGGATCGGCTGCTAAAGCTA" ) );
	mm.build( X );
	string_t s = "AC";
	std::cout << "p(" << s << ")=" << mm.emittedProbability( s ) << "\n";

	return 0;
}
