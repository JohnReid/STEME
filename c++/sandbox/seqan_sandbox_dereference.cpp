/**
 * Copyright John Reid 2009
 */

#include <seqan/index.h>

using namespace seqan;

typedef String< Dna5 >                                        string_t;              /**< A string of Dna5. */
typedef Iterator< string_t >::Type                            string_it;             /**< An iterator over our strings. */
typedef Value< string_t >::Type                               string_value_t;        /**< The value type of our string. */
typedef Reference< string_t >::Type                           string_ref_t;          /**< The reference type of our string. */

int
main( int argc, char * argv[] ) {
	//
	// Test for assigning dereferenced iterator
	//
	string_t                            seq( "acgt" );
	string_it                           seq_start = begin( seq );
	string_value_t                      c = *seq_start; // this works fine
	string_ref_t                        seq_ref = *seq_start; // this works fine
	unsigned long                       d = *seq_start; //!!!!!!!!! causes seg fault

	return 0;
}
