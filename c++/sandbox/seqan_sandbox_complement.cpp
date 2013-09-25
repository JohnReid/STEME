/**
 * Copyright John Reid 2009
 */

#include <seqan/index.h>

int
main( int argc, char * argv[] )
{
	using namespace seqan;

	//
	// Test for FunctorComplement< Dna5 >
	//
	FunctorComplement< Dna5 > comp;
	std::string bases = "acgtn";
	for( size_t i = 0; 5 != i; ++i ) {
		Dna5 b = bases[ i ];
		std::cout << b  << " : " << comp( b ) << " : " << comp( comp( b ) ) << "\n";
	}

	return 0;
}
