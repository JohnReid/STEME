/**
 * Copyright John Reid 2010
 */

#include <seqan/sequence.h>
#include <seqan/modifier.h>
#include <seqan/file.h>
#include <iostream>

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

/// Print a string and its reverse complement
template< typename T >
void
show_rev_comp( T t ) {
	typedef typename rev_comp< T >::type rev_comp_t;

	std::cout << "Original:           " << t << "\n";
	std::cout << "Reverse complement: " << rev_comp_t( t ) << "\n";
}

int
main( int argc, char * argv[] )
{
	using namespace seqan;

	typedef String< Dna5 >                    string_t;
	typedef Infix< string_t >::Type           infix_t;
	typedef Prefix< string_t >::Type          prefix_t;
	typedef Suffix< string_t >::Type          suffix_t;

	string_t string = "ACGTNACGTN";
	show_rev_comp( string );

	infix_t infix = infixWithLength( string, 3, 5 );
	show_rev_comp( infix );

	prefix_t pre = prefix( string, 5 );
	//show_rev_comp( pre ); // <-- CAUSES COMPILATION ERROR

	suffix_t suf = suffix( string, 5 );
	//show_rev_comp( suf ); // <-- CAUSES COMPILATION ERROR

	return 0;
}
