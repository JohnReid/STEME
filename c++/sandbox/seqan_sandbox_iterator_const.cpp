/**
 * Copyright John Reid 2010
 */

#include <seqan/sequence.h>
#include <seqan/modifier.h>

/// Defines a type to represent the reverse complement of the template string type.
template< typename StringT >
struct rev_comp {

	typedef typename seqan::Value< StringT >::Type   string_char; ///< The type of a character in the string.

	/// The type of the reverse complement of the string.
	typedef seqan::ModifiedString<
		seqan::ModifiedString<
			StringT,
			seqan::ModView< seqan::FunctorComplement< string_char > >
		>,
		seqan::ModReverse
	> type;
};

int
main( int argc, char * argv[] ) {
	using namespace seqan;

	typedef String< Dna5 >                                string_t;
	typedef rev_comp< string_t >::type                    rev_comp_t;

	string_t                                        str( "ACGTNACGTNACGTN" );
	rev_comp_t                                      rev_comp( str );
	begin( str );
	begin( rev_comp );

	const string_t                                  const_str( "ACGTNACGTNACGTN" );
	const rev_comp_t                                const_rev_comp( str );
	begin( const_str );
	//begin( const_rev_comp ); // <-- CAUSES COMPILATION ERROR
}

