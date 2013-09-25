/**
 * Copyright John Reid 2011
 */

#include <boost/progress.hpp>
#include <iostream>
#include <cmath>
#include <vector>

#define _STEM_LOG_BITS 64.0

/// Log(x) for some x close to 0.
#define _STEM_LOG_NEAR_ZERO -1e100

/// log(x + y) where log(x) >= log(y).
#define _STEM_LOG_X_PLUS_Y(logx,logy) ( \
	( (logy) <= _STEM_LOG_NEAR_ZERO || ( (logx) - (logy) ) > _STEM_LOG_BITS ) \
		? (logx) \
		: (logx) + log1p( std::exp( (logy) - (logx) ) ) \
)

/// log(x + y) where log(x) and log(y) are available.
#define STEM_LOG_X_PLUS_Y(logx,logy) ( \
	( (logx) > (logy) ) \
		? _STEM_LOG_X_PLUS_Y( (logx), (logy) ) \
		: _STEM_LOG_X_PLUS_Y( (logy), (logx) ) \
)

int
main( int argc, char * argv[] ) {

	const double log_a = -.3;
	const double log_b = -.4;

	const size_t num_cycles = 1e10;
	std::vector< double > results( 1000 );

	{
		boost::progress_timer t( std::cout );  // start timing
		for( size_t i = 0; num_cycles != i; ++i ) {
			results[ i % 1000 ] = std::exp( log_a - STEM_LOG_X_PLUS_Y( log_a, log_b ) );
		}
	}

	{
		boost::progress_timer t( std::cout );  // start timing
		for( size_t i = 0; num_cycles != i; ++i ) {
			results[ i % 1000 ] = std::exp( log_a ) / ( std::exp( log_a ) + std::exp( log_b ) );
		}
	}

	return 0;
}
