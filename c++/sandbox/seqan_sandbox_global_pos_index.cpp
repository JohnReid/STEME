/**
 * Copyright John Reid 2011
 *
 * \file
 * Code to check if using enhanced suffix arrays with SAValue's of global positions works
 */

#include <seqan/index.h>
#include <seqan/sequence.h>

using namespace seqan;

struct global {};

// specialise for our index
namespace seqan {

template< typename TText >
struct SAValue< Index< TText, IndexEsa< global > > >
{
	typedef unsigned Type;
};

} // namespace seqan



typedef String< Dna5 >                            string_t;              ///< A string of Dna5.
typedef StringSet< string_t >                     string_set_t;          ///< StringSet type.
typedef Index< string_set_t, IndexEsa< global > > index_t;

int
main( int argc, char * argv[] ) {

	string_set_t sequences;
	index_t index( sequences );
	//Iterator< index_t, TopDown<> >::Type it( index ); // <<<<<<<<<<<<<<<<<< this line causes compilation failure

	return 0;
}
