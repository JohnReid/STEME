/**
 * Copyright John Reid 2011, 2016
 *
 * \file Code to calculate p-values.
 */


#ifndef _PVALUES_JR_22AUG2011_DEFS_H_
#define _PVALUES_JR_22AUG2011_DEFS_H_

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/range.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/assert.hpp>
#include <boost/multi_array.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/round.hpp>

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iostream>
#include <sstream>


// turn off pvalue checks for NDEBUG builds
#ifdef NDEBUG
# define PVALUE_CHECKS_DISABLED
#endif


// define PVALUE_ASSERT as BOOST_ASSERT if not disabled
#ifdef PVALUE_CHECKS_DISABLED
# define PVALUE_ASSERT( x )
#else // PVALUE_CHECKS_DISABLED
# define PVALUE_ASSERT( x ) BOOST_ASSERT( x )
#endif //PVALUE_CHECKS_DISABLED


/** Creates a string representation of the streamed message. */
#define PVALUE_MAKE_STRING( msg ) \
    ( ( std::ostringstream & )( std::ostringstream() << std::dec << msg ) ).str()


/**
 * Code to calculate p-values.
 */
namespace pvalues {


// forwards
template< typename Spec > struct log_arithmetic;


/** Namespace for specifications. */
namespace spec {

/** Default specification tag. */
struct Default { };

/** Float meta-function. */
template< typename Spec >
struct Float {
    typedef double type;
};

/** Log arithmetic meta-function. */
template< typename Spec >
struct LogArithmetic {
    typedef log_arithmetic< Spec > type;
};

/** Log on/off meta-function. */
template< typename Spec >
struct LogOn {
    static const bool value = false;
};

} // namespace spec



/// Implements log lookup operations
template< typename T >
struct log_lookup {

    typedef std::vector< T > vec;
    typedef typename vec::size_type size_t;

    static BOOST_CONSTEXPR_OR_CONST int D_LOGADD = 30; // y/x < 1e-D discard
    static BOOST_CONSTEXPR_OR_CONST int G_LOGADD = 100; // step function look-up every 1/G
    static BOOST_CONSTEXPR_OR_CONST int D_LOGSUB = 17; // above 16.25.. the func is zero
    static BOOST_CONSTEXPR_OR_CONST int G_LOGSUB = 200;
    static BOOST_CONSTEXPR_OR_CONST double SUB_INTERP_THRESHOLD = .02; // too steep to interpolate on

    static inline T smallest() { return - std::numeric_limits< T >(); }

    size_t d_logadd, d_logsub;
    vec logaddf, logsubf;

    log_lookup()
    : d_logadd( D_LOGADD * std::log( 10. ) * G_LOGADD )
    , d_logsub( D_LOGSUB * std::log( 10. ) * G_LOGSUB )
    , logaddf( d_logadd + 1 )
    , logsubf( d_logsub + 1 )
    {
        // initialise logaddf
        for( size_t i = 0; d_logadd + 1 != i; ++i ) {
            logaddf[ i ] = boost::math::log1p( std::exp( -T( i ) / G_LOGADD ) );
        }

        // initialise logsubf
        logsubf[ 0 ] = -std::numeric_limits< T >::max();
        for( size_t i = 1; d_logsub + 1 != i; ++i ) {
            logsubf[ i ] = boost::math::log1p( -std::exp( -T( i ) / G_LOGSUB ) );
        }
    }

    /// \return log(x+y) when given log(x) and log(y)
    inline
    T
    log_add( T log_x, T log_y ) {
        if( log_x == smallest() )
            return log_y;
        else if( log_y == smallest() ) {
            return log_x;
        } else {
            // make log_x bigger than log_y by swapping if necessary
            if( log_x < log_y ) {
                std::swap( log_x, log_y );
            }
            // get lookup index
            const size_t i = size_t( ( log_x - log_y ) * G_LOGADD );
            if( i < d_logadd ) { // linear interpolation
                log_x +=
                    ( i + 1 - ( log_x - log_y ) * G_LOGADD ) * logaddf[ i ]
                    + ( ( log_x - log_y ) * G_LOGADD - i ) * logaddf[ i + 1 ]
                    ;
            }
        }
    }

    /// \return log(x-y) when given log(x) and log(y)
    inline
    T
    log_sub( T log_x, T log_y ) {
        const T diff = log_x - log_y;
        if( diff <= 0. ) {
            // underflow++; // we're not using underflow counter
            return smallest();
        } else if( log_y == smallest() ) {
            return log_x;
        } else {
            if( diff < SUB_INTERP_THRESHOLD ) { // too steep to interpolate on
                return log_x + boost::math::log1p( -std::exp( -diff ) );
            } else {
                // lookup index
                const size_t i = size_t( ( diff ) * G_LOGSUB );
                if( i < d_logsub ) // linear interpolation
                    return
                        log_x
                        + ( i + 1 - diff * G_LOGSUB ) * logsubf[ i ]
                        + ( diff * G_LOGSUB - i ) * logsubf[ i + 1 ]
                        ;
            }
        }
    }
};



/**
 * \return log(x+y) when given log(x) and log(y).
 *
 * Makes one std::exp() call and one boost::math::log1p() call
 */
inline
double
log_add( double log_x, double log_y ) {
    // make log_x the larger of the two.
    if( log_x < log_y ) {
        std::swap( log_x, log_y );
    }
    return log_x + boost::math::log1p( std::exp( log_y - log_x ) );
}


/**
 * Defines a system to perform arithmetic operations on logarithms.
 *
 * Can be specialised.
 */
template< typename Spec = spec::Default >
struct log_arithmetic {
    typedef typename spec::Float< Spec >::type float_t; ///< Floating point type.
    typedef log_arithmetic< Spec > self_t; ///< This type.
    typedef float_t log_t; ///< Type that stores logarithms.
    typedef std::vector< log_t > log_vec; ///< Vector of logarithms.
    typedef float_t real_t; ///< Type that stores real values.
    typedef unsigned count_t; ///< Type that stores counts.

protected:
    log_vec _log_n; ///< Logs of counts.
    log_vec _log_fact_n; ///< Logs of factorials.

public:
    /// \return representation of log(0).
    inline
    log_t
    log_zero() {
        return log( 0u );
    }

    /// \return log(x).
    inline
    log_t
    log( real_t x ) {
        return std::log( x );
    }

    /// \return exp(log(x)) = x
    inline
    real_t
    exp( log_t log_x ) {
        return std::exp( log_x );
    }

    /// \return log(n).
    inline
    log_t
    log( count_t n ) {
        PVALUE_ASSERT( n < _log_n.size() );
        return _log_n[ n ];
    }

    /// \return log(n!).
    inline
    log_t
    log_fact( count_t n ) {
        PVALUE_ASSERT( n < _log_fact_n.size() );
        return _log_fact_n[ n ];
    }

    /// \return log(x+y) given log(x) and log(y).
    inline
    log_t
    log_add( log_t log_x, log_t log_y ) {
        const log_t max_x_y = std::max( log_x, log_y );
        if( boost::math::isinf( log_x ) || boost::math::isinf( log_x ) ) {
            return max_x_y;
        }
        return max_x_y + boost::math::log1p( std::exp( - std::fabs( log_x - log_y ) ) );
    }

    /// \return log(x-y) given log(x) and log(y). Pre-condition: log(x) > log(y).
    inline
    log_t
    log_subtract( log_t log_x, log_t log_y ) {
        PVALUE_ASSERT( log_x > log_y );
        if( boost::math::isinf( log_y ) ) {
            return log_x;
        }
        return log_x + boost::math::log1p( - std::exp( log_y - log_x ) );
    }

    /// \return Whether the log value is finite
    inline
    bool
    log_isfinite( log_t log_x ) {
        return boost::math::isfinite( log_x );
    }

    /// \return Maximum n
    inline
    size_t
    max_n() const {
        return _log_n.size() - 1;
    }

    /// Constructor.
    log_arithmetic( count_t max_n ) {
        using namespace boost::lambda;

        // calculate log(n_i) for n_i in [0, max_n]
        _log_n.reserve( max_n + 1 );
        _log_n.push_back( - std::numeric_limits< log_t >::max() ); // set log(0) to smallest value.
        typename self_t::log_t ( self_t::* log_fn_ptr )( typename self_t::real_t ) = &self_t::log;
        std::copy(
            boost::make_transform_iterator( boost::counting_iterator< count_t >( count_t( 1 ) ), bind( log_fn_ptr, this, boost::lambda::_1 ) ),
            boost::make_transform_iterator( boost::counting_iterator< count_t >( max_n + 1    ), bind( log_fn_ptr, this, boost::lambda::_1 ) ),
            std::back_inserter( _log_n )
        );

        // calculate log(n_i!)
        _log_fact_n.reserve( _log_n.size() );
        _log_fact_n.push_back( 0. ); // set log(0!) to 0.
        for( count_t i = 1; _log_n.size() != i; ++i ) {
            _log_fact_n.push_back( std::accumulate( _log_n.begin() + 1, _log_n.begin() + i + 1, 0. ) );
        }
    }
};




/// Does binary search. Finds first index, i, where Fn(i) is true.
template< typename Fn >
struct binary_search {

    /// Construct.
    binary_search( Fn & fn ) : _fn( fn ) { }

    /// Perform the binary search.
    template< typename IndexT >
    inline
    IndexT
    operator()( IndexT min, IndexT max ) {
        // check extreme values
        PVALUE_ASSERT( min <= max );
        IndexT result;
        if( max == min ) {
            result = max;
        } else if( ! _fn( max - 1 ) ) {
            result = max;
        } else if( _fn( min ) ){
            result = min;
        } else {
            // do binary search
            PVALUE_ASSERT( ! _fn( min ) );
            PVALUE_ASSERT( _fn( max - 1 ) );
            result = unchecked_search( min, max );
        }

        // check that we found one of the extrema or that our functional does go
        // from false to true at result
        PVALUE_ASSERT( result == min || ! _fn( result - 1 ) );
        PVALUE_ASSERT( result == max || _fn( result ) );

        return result;
    }

protected:
    Fn & _fn;

    // binary search without checking values of end-points
    template< typename IndexT >
    inline
    IndexT
    unchecked_search( IndexT min, IndexT max ) {
        const IndexT diff = max - min;
        if( 1 == diff ) {
            return max;
        } else {
            PVALUE_ASSERT( ! _fn( min ) );
            const IndexT mid_point = min + diff / 2;
            if( _fn( mid_point ) ) {
                return unchecked_search( min, mid_point );
            } else {
                return unchecked_search( mid_point, max );
            }
        }
    }
};




/**
 * Do a binary search.
 */
template<
    typename Fn,
    typename IndexT
>
IndexT
do_binary_search( Fn fn, IndexT min, IndexT max ) {
    return binary_search< Fn >( fn )( min, max );
}



/**
 * Negate a binary function.
 */
template< typename Fn, typename ArgT  >
struct unary_negate
: std::unary_function< ArgT, bool >
{
    Fn & _fn;

    unary_negate( Fn & fn ) : _fn( fn ) { }

    bool
    operator()( ArgT x ) {
        return ! _fn( x );
    }
};


template<
    typename ArgT,
    typename Fn
>
unary_negate< Fn, ArgT >
make_unary_negate( Fn & fn ) {
    return unary_negate< Fn, ArgT >( fn );
}



/// Normalise the log cumulative distribution of the log probability mass function.
template< typename CmfRange >
typename boost::range_value< CmfRange >::type
normalise_log_cumulative( CmfRange & log_cmf ) {
    typedef typename boost::range_value< CmfRange >::type float_t;
    const float_t offset = - *boost::begin( log_cmf );
    BOOST_FOREACH( float_t & lc, log_cmf ) {
        lc += offset;
    }
    return offset;
}


/// Calculate the log cumulative distribution of the log probability mass function.
template< typename PmfRange >
boost::shared_array< typename boost::range_value< PmfRange >::type >
calculate_log_cumulative( const PmfRange & log_pmf, bool normalise = true ) {
    typedef typename boost::range_value< PmfRange >::type float_t;
    const size_t Q = boost::end( log_pmf ) - boost::begin( log_pmf );
    boost::shared_array< double > result( new float_t[ Q ] );
    float_t log_cdf = std::log( 0. );
    for( int q = Q; 0 != q; ) {
        log_cdf = log_add( log_cdf, log_pmf[ --q ] );
        result[ q ] = log_cdf;
    }
    if( normalise ) {
        boost::iterator_range< double * > range( result.get(), result.get() + Q );
        const double norm_adjustment = normalise_log_cumulative( range );
        // std::cout << "Normalisation adjustment " << norm_adjustment << "\n";
    }
    return result;
}







} // namespace pvalues

#endif // _PVALUES_JR_22AUG2011_DEFS_H_
