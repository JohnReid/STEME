/**
 * Copyright John Reid 2011
 *
 * \file Code to parse p-value test cases */


#include "pvalue_test_defs.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>


using namespace boost;
using namespace std;

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;


/**
 * Our test parser.
 */
template< typename Iterator >
optional< p_value_test_case >
parse_IC_line( Iterator first, Iterator last )
{
	using qi::double_;
	using qi::_1;
    using qi::uint_;
    using qi::double_;
    using qi::phrase_parse;
    using ascii::space;
    using boost::phoenix::ref;

    optional< p_value_test_case > result;

    /// Will look like "N=2; IC=1.57184,1.57184,1.57184,1.22985,1.22985,1.57184,1.57184,1.22985"
    size_t N;
    if( phrase_parse( first, last, "N=" >> uint_ >> "; IC=", space, N ) ) {
    	result.reset( p_value_test_case( N ) );
        if( ! phrase_parse( first, last, double_ % ",", space, result->ICs ) ) {
        	result = optional< p_value_test_case >(); // reset
        }
    }

    return result;
}


void
read_p_value_test_cases(
    ifstream & in,
    p_value_test_case_vec & arguments,
    size_t min_N,
    size_t max_N,
    size_t max_cases
) {
	string line;
	while( ! in.eof() ) {
		getline( in, line );
		if( in.eof() || ! in ) break;
		optional< p_value_test_case > args = parse_IC_line( line.begin(), line.end() );
		if( args ) {
            if( ! min_N || args->N >= min_N ) {
                if( ! max_N || args->N <= max_N ) {
                    arguments.push_back( *args );
                    if( max_cases && arguments.size() >= max_cases ) {
                        break;
                    }
                }
			}
		}
	}
}

