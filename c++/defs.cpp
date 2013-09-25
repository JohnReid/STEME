/** Copyright John Reid 2008
*/




#include "steme_defs_pch.h"
#include <steme/defs.h>

//I/O for matrix types
namespace std {

/** Stream a matrix row. */
std::ostream &
operator<<( std::ostream & os, steme::double_matrix_t::const_reference row )
{
	std::copy( row.begin(), row.end(), std::ostream_iterator< steme::double_matrix_t::const_reference::value_type >( os, "\t" ) );
	return os;
}



/** Stream a matrix. */
std::ostream &
operator<<( std::ostream & os, const steme::double_matrix_t & m )
{
	//os << m[0] << "\n";
	std::copy( m.begin(), m.end(), std::ostream_iterator< steme::double_matrix_t::const_reference >( os, "\n" ) );
	return os;
}






} //namespace std







namespace steme {


//precalculate some constants
const double log_quarter = std::log( 0.25 ); /**< log(.25) */
const double log_half = std::log( 0.5 );     /**< log(.5) */
const double log_95 = std::log( 0.95 );      /**< log(.95) */
const double log_05 = std::log( 0.05 );      /**< log(.05) */




double p_true_from_log_odds( double log_odds )
{
	const double odds = std::exp( log_odds );
	return odds / ( odds + 1. );
}

double log_odds_from_p_true( double p_true )
{
	return std::log( p_true ) - std::log( 1. - p_true );
}

double_matrix_ptr reverse_complement_of_matrix( const double_matrix_t & m ) {
	const unsigned K = m.shape()[0];
	const unsigned B = m.shape()[1]; //normally 4
	double_matrix_ptr result( new double_matrix_t( boost::extents[K][B] ) );
	for( unsigned i = 0; K != i; ++i ) {
		for( unsigned j = 0; B != j; ++j ) {
			(*result)[K-1-i][B-1-j] = m[i][j];
		}
	}
	return result;
}



void normalise_counts( double_matrix_t & counts )
{
	BOOST_FOREACH( double_matrix_t::reference r, counts ) {
		double factor = 1.0 / std::accumulate( r.begin(), r.end(), 0.0 );
		BOOST_FOREACH( double_matrix_t::reference::reference e, r ) {
			e *= factor;
		}
	}
}



double_matrix_ptr get_log_likelihoods_from_nucleo_dist( const double_matrix_t & dist )
{
	double_matrix_ptr result( new double_matrix_t( boost::extents[dist.shape()[0]][dist.shape()[1]] ) );
	for( unsigned i = 0; dist.shape()[0] != i; ++i ) {
		for( unsigned j = 0; dist.shape()[1] != j; ++j ) {
			(*result)[i][j] = std::log( dist[i][j] );
		}
	}
	return result;
}

} //namespace steme
