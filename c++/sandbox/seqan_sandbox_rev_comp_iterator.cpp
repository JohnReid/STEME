/**
 * Copyright John Reid 2010
 */

#include <seqan/sequence.h>
#include <seqan/modifier.h>
#include <seqan/file.h>
#include <iostream>


int main( int argc, char * argv [] ) {
	using namespace seqan;

	typedef Dna5                                                            dna_t;
	typedef String< dna_t >                                                 string_t;
	typedef Iterator< string_t >::Type                                      iterator_t;
	typedef ModifiedIterator< iterator_t, ModReverse >                      rev_iterator_t;
	typedef ModifiedIterator< rev_iterator_t, FunctorComplement< dna_t > >  rev_comp_iterator_t;

	string_t str = "ACGTNAAGGTTCC";
	::std::cout << "Original string: " << str << "\n";

	rev_iterator_t rev_begin( end( str ) );
	rev_iterator_t rev_end( begin( str ) );
	::std::cout << "Reverse string:  " << str;
	while( rev_begin != rev_end ) {
		::std::cout << *rev_begin;
		++rev_begin;
	}
	::std::cout << "\n";


//	rev_comp_iterator_t rev_comp_begin( rev_iterator_t( end( str ) ) );
//	rev_comp_iterator_t rev_comp_end( rev_iterator_t( begin( str ) ) );
//
//	::std::cout << "Rev comp string: " << str;
//	while( rev_comp_begin != rev_comp_end ) {
//		::std::cout << *rev_comp_begin;
//		++rev_comp_begin;
//	}
//	::std::cout << "\n";

	return 0;
}
