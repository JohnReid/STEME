/** Copyright John Reid 2011, 2012, 2013, 2014
 *
 * \file
 * \brief Basic definitions and includes for STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_DEFS_H_
#define STEME_JR_13AUG2011_DEFS_H_

#include <myrrh/assert.h>
#include <myrrh/seqan_boost_range_adaptors.h>

#include <seqan/index.h>
#include <seqan/version.h>

#include <boost/foreach.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/multi_array.hpp>
#include <boost/operators.hpp>
#include <boost/optional.hpp>
#include <boost/range/iterator_range.hpp>


#include <algorithm>
#include <vector>
#include <sstream>
#include <string>


#if SEQAN_VERSION_MAJOR > 1 || SEQAN_VERSION_MINOR > 3
# define SEQAN_PACK Pack
# define SEQAN_BITPACKED BitPacked
#else
# define SEQAN_PACK Compressed
# define SEQAN_BITPACKED BitCompressed
#endif

/// log(.5)
#define STEME_LOG_HALF -0.69314718055994529

#define _STEME_LOG_BITS 64.0

/// Log(x) for some x close to 0.
#define _STEME_LOG_NEAR_ZERO -1e100

/// log(x + y) where log(x) >= log(y).
#define _STEME_LOG_X_PLUS_Y(logx,logy) ( \
	( (logy) <= _STEME_LOG_NEAR_ZERO || ( (logx) - (logy) ) > _STEME_LOG_BITS ) \
		? (logx) \
		: (logx) + log1p( std::exp( (logy) - (logx) ) ) \
)

/// log(x + y) where log(x) and log(y) are available.
#define STEME_LOG_X_PLUS_Y(logx,logy) ( \
	( (logx) > (logy) ) \
		? _STEME_LOG_X_PLUS_Y( (logx), (logy) ) \
		: _STEME_LOG_X_PLUS_Y( (logy), (logx) ) \
)




#if defined( MS_HAS_VALGRIND ) && ! defined( NDEBUG )
# define MS_USING_VALGRIND 1
# include <valgrind/memcheck.h>
#else //defined( HAS_VALGRIND ) && ! defined( NDEBUG )
# define MS_USING_VALGRIND 0 ///< Whether we are using valgrind or not.
#endif //defined( HAS_VALGRIND ) && ! defined( NDEBUG )




#if MS_USING_VALGRIND
# define MS_ASSERT_VALUE_IS_DEFINED( value ) \
	if( VALGRIND_CHECK_VALUE_IS_DEFINED( value ) ) { \
		MYRRH_ASSERT( ! "Valgrind: " #value " is not defined." ); \
	}
#else //MS_USING_VALGRIND
# define MS_ASSERT_VALUE_IS_DEFINED( value ) ///< Assert that value is defined (uses VALGRIND)
#endif //MS_USING_VALGRIND



#if !defined NDEBUG && !defined MYRRH_DISABLE_ASSERTS

#define MS_DO_DEBUG_CHECKS

/** Check that isnan(x) is false. */
#define MS_CHECK_IS_NUMBER( x, name ) \
	if( ::boost::math::isnan( x ) ) { \
		throw std::logic_error( name " is NAN" ); \
	}

/** Check that isfinite(x) is true. */
#define MS_CHECK_IS_FINITE( x, name ) \
	if( ! ::boost::math::isfinite( x ) ) { \
		throw std::logic_error( name " is not finite" ); \
	}

/** Check that x is inside given bounds is true. */
#define MS_CHECK_BOUNDS( x, name, lower, upper ) \
	if( ( x ) < lower || upper < ( x ) ) { \
		throw std::logic_error( "Invalid value of " name " outside bounds" ); \
	}

#else //!defined NDEBUG && !defined MYRRH_DISABLE_ASSERTS

#undef MS_DO_DEBUG_CHECKS

#define MS_CHECK_IS_NUMBER( x, name )
#define MS_CHECK_IS_FINITE( x, name )
#define MS_CHECK_BOUNDS( x, name, lower, upper )

#endif //!defined NDEBUG && !defined MYRRH_DISABLE_ASSERTS

/** Creates a string representation of the streamed message. */
#define MS_MAKE_STRING( msg ) \
	( ( std::ostringstream & )( std::ostringstream() << std::dec << msg ) ).str()


//#define MS_NO_DEBUG_STRINGS // define this to disable MS_DEBUG_STRING macros even in debug builds
#if defined(NDEBUG) || defined(MS_NO_DEBUG_STRINGS)
/**
In debug builds define a std::string variable with name, var_name and copy seqan::String, s, into it.
This is useful in order to inspect the values of seqan strings in the debugger.
*/
# define MS_DEBUG_STRING_VAR( var_name, s )
/**
In debug builds define a std::string variable with name, var_name and copy the range into it.
This is useful in order to inspect the values of seqan strings in the debugger.
*/
# define MS_DEBUG_STRING_RANGE( var_name, range )

#else //defined(NDEBUG) || defined(MS_NO_DEBUG_STRINGS)

/**
In debug builds define a std::string variable with name, var_name and copy seqan::String, s, into it.
This is useful in order to inspect the values of seqan strings in the debugger.
*/
# define MS_DEBUG_STRING_VAR( var_name, s ) std::string var_name; seqan::append( var_name, s )
/**
In debug builds define a std::string variable with name, var_name and copy the range into it.
This is useful in order to inspect the values of seqan strings in the debugger.
*/
# define MS_DEBUG_STRING_RANGE( var_name, range ) \
	std::string var_name; \
	std::copy( boost::begin( range ), boost::end( range ), std::back_insert_iterator< std::string >( var_name ) )

#endif //defined(NDEBUG) || defined(MS_NO_DEBUG_STRINGS)



/**
 * Macro from http://dbp-consulting.com/tutorials/SuppressingGCCWarnings.html to
 * control GCC warnings
 */
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402
# define GCC_DIAG_STR(s) #s
# define GCC_DIAG_JOINSTR(x,y) GCC_DIAG_STR(x ## y)
# define GCC_DIAG_DO_PRAGMA(x) _Pragma (#x)
# define GCC_DIAG_PRAGMA(x) GCC_DIAG_DO_PRAGMA(GCC diagnostic x)
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#  define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(push) \
    GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x))
#  define GCC_DIAG_ON(x) GCC_DIAG_PRAGMA(pop)
# else
#  define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x))
#  define GCC_DIAG_ON(x)  GCC_DIAG_PRAGMA(warning GCC_DIAG_JOINSTR(-W,x))
# endif
#else
# define GCC_DIAG_OFF(x)
# define GCC_DIAG_ON(x)
#endif


/**
 * It seems that gcc 4.6 and 4.7 name this warning differently
*/
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 407
# define GCC_UNINIT_WARNING maybe-uninitialized
#else
# define GCC_UNINIT_WARNING uninitialized
#endif



namespace steme {

typedef std::vector< double >             double_vec;                ///< A vector of doubles.
typedef std::vector< double_vec >         double_vec_vec;            ///< A vector of vector of doubles.
typedef boost::multi_array< double, 1 >   double_array;              ///< An array of doubles.

// precalculated constants
extern const double log_quarter;
extern const double log_half;
extern const double log_95;
extern const double log_05;



/**
 * Convert a likelihood ratio to a probability.
 */
template< typename FPT >
FPT
lr_to_prob( FPT lr ) {
    return lr / ( lr + FPT( 1 ) );
}


/**
 * Convert a likelihood ratio to a probability.
 */
template< typename FPT >
FPT
prob_to_lr( FPT p ) {
    return p / ( FPT( 1 ) - p );
}


/**
 * Do two W-mers overlap?
 */
inline
bool
wmer_overlap( int pos1, int pos2, int W ) {
    return std::abs( pos1 - pos2 ) < W;
}


/**
 * Shrinks a vector to empty to release the memory it is holding.
 */
template< typename T, class Alloc >
void
shrink_to_empty( std::vector< T, Alloc > & v )
{
    // swap with empty vector
    std::vector< T, Alloc >( v.begin(), v.end(), v.get_allocator() ).swap( v );
}


/**
 * Return a std::string version of a seqan string.
 */
template< typename SeqanStrT >
std::string
std_string_from_seqan( SeqanStrT seqan_str ) {
	std::string result;
	seqan::assign( result, seqan_str );
	return result;
}


/** Check the exponentiated log probabilities sum close to 1. */
template< typename LogProbRange >
bool
log_probs_sum_to_1( const LogProbRange & log_probs, double tolerance = .01 ) {
	typedef typename boost::range_value< LogProbRange >::type value_t;
	value_t sum = 0.;
	BOOST_FOREACH( value_t p, log_probs ) {
		sum += std::exp( p );
	}
	return (1. - tolerance) < sum && sum < (1. + tolerance);
}

/** Check the probabilities sum close to 1. */
template< typename ProbRange >
bool
probs_sum_to_1( const ProbRange & probs, double tolerance = .01 ) {
	typedef typename boost::range_value< ProbRange >::type value_t;
	value_t sum = 0.;
	BOOST_FOREACH( value_t p, probs ) {
		sum += p;
	}
	return (1. - tolerance) < sum && sum < (1. + tolerance);
}



/// Descend a suffix tree/array index to a given depth, W
template< typename Derived >
struct descend_to_depth {

	size_t W; ///< Depth to descend to.

	/// Constructor.
	descend_to_depth( size_t W ) : W( W ) { }

	/// Descend the index to the requested depth.
	template< typename Index >
	void
	descend_index( Index & index ) {
		typedef typename seqan::Iterator< Index, seqan::TopDown<> >::Type top_down_it;
		visit_node( top_down_it( index ) );
	}

protected:
	template< typename It >
	void
	visit_node( It it ) {

		using namespace seqan;

		const size_t w = std::min( repLength( it ), W );

		MS_DEBUG_STRING_VAR( w_mer, prefix( representative( it ), w ) );

		if( W == w ) {

			// we have reached desired depth so perform action
			static_cast< Derived * >( this )->handle_node( it );

		} else {

			// continue descending the suffix tree, so go down and across
			if( goDown( it ) ) {
				visit_node( it );
				while( goRight( it ) ) {
					visit_node( it );
				}
			}

		}
	}
};


/** Normally used for background frequencies. Enforces complement constraint, i.e. p(A)=p(T) and p(C)=p(G). */
struct zero_order_frequencies {

	typedef boost::multi_array< double, 1 >              double_array;       ///< An array of doubles.
	typedef boost::shared_ptr< zero_order_frequencies >  ptr;                ///< A pointer to a zero_order_frequencies object.

	double_array               dist;          ///< Background frequency distribution.
	double_array               dist_logs;     ///< Logs of background frequency distribution.
	size_t                     total_counts;  ///< Holds a record of the total counts used to make the distribution.

private:
	/// Default constructor.
	zero_order_frequencies()
	: dist( boost::extents[ 4 ] )
	, dist_logs( boost::extents[ 4 ] )
	{ }

	/// Normalise
	void
	normalise() {
		total_counts = std::accumulate( dist.begin(), dist.end(), double( 0 ) );
		const double p_a_or_t = ( dist[0] + dist[3] ) / total_counts;
		dist[0] = dist[3] = p_a_or_t / 2.;
		dist[1] = dist[2] = (1. - p_a_or_t) / 2.;
		MYRRH_ASSERT( probs_sum_to_1( dist ) );
		calculate_logs();
	}


	/// calculate logs of background dist
	void calculate_logs() {
		using boost::lambda::bind;
		double( *log_fn )( double ) = std::log;
		std::transform( dist.begin(), dist.end(), dist_logs.begin(), bind( log_fn, boost::lambda::_1 ) );
		MYRRH_ASSERT( log_probs_sum_to_1( dist_logs ) );
	}


public:
    /// Constructor.
    template< typename It >
    zero_order_frequencies( It begin, It end )
    : dist( boost::extents[ 4 ] )
    , dist_logs( boost::extents[ 4 ] )
    {
        std::copy(
            begin,
            end,
            dist.begin()
        );
        normalise();
    }


	/// Return another zero_order_frequencies object with the pseudo-counts added.
	zero_order_frequencies
	add_pseudo_counts( double pseudo_counts ) const {
		zero_order_frequencies result;
		result.total_counts = total_counts + 4. * pseudo_counts;
		result.dist[0] = result.dist[3] = (this->dist[0] * this->total_counts + pseudo_counts) / result.total_counts;
		result.dist[1] = result.dist[2] = (this->dist[1] * this->total_counts + pseudo_counts) / result.total_counts;
		MYRRH_ASSERT( probs_sum_to_1( result.dist ) );
		result.calculate_logs();
		return result;
	}
};



/// Holds a pair of optionals of the templated type. Implements comparison operators.
template< typename T >
struct optional_pair
: boost::less_than_comparable< optional_pair< T > >
, boost::equality_comparable< optional_pair< T > >
{
	typedef boost::optional< T >    optional_t;      ///< Optional type.

	optional_t                      first;           ///< First optional.
	optional_t                      second;          ///< Second optional.

	/** Default constructor. */
	optional_pair() {
	}

	/** Constructor. */
	optional_pair( optional_t first, optional_t second )
	: first( first )
	, second( second )
	{
	}

	/// True if either optional has a value.
	operator bool() const {
		return first || second;
	}

	/// The greater value.
	optional_t
	greater() const {
		if( optional_less_than( first, second ) ) return second;
		return first;
	}

	/// The lesser value.
	optional_t
	lesser() const {
		if( optional_less_than( first, second ) ) return first;
		return second;
	}

	/// Comparison.
	bool
	operator<( optional_pair< T > rhs ) const {
		if( this->greater() < rhs.greater() ) return true;
		if( rhs.greater() < this->greater() ) return false;
		return this->lesser() < rhs.lesser();
	}

	/// Equality.
	bool
	operator==( optional_pair< T > rhs ) const {
		return this->first == rhs.first && this->second == rhs.second;
	}

	/// Equality for optional type
	static
	bool
	optional_equality( optional_t lhs, optional_t rhs ) {
		if( lhs ) {
			if( rhs ) return *lhs == *rhs;
			return false;
		} else {
			if( rhs ) return false;
			return true;
		}
	}

	/// Comparison for optional type
	static
	bool
	optional_less_than( optional_t lhs, optional_t rhs ) {
		if( lhs ) {
			if( rhs ) return *lhs < *rhs;
			return true;
		} else {
			if( rhs ) return true;
			return false;
		}
	}
};


/// Scale the components of the optional pair.
template< typename T >
optional_pair< T > &
scale_optional_pair( optional_pair< T > & pair, T scale)
{
	if( pair.first ) *pair.first *= scale;
	if( pair.second ) *pair.second *= scale;
	return pair;
}


/// Add the 2 optional pairs.
template< typename T >
optional_pair< T > &
operator+=( optional_pair< T > & lhs, const optional_pair< T > & rhs )
{
	if( lhs.first || rhs.first ) lhs.first = lhs.first.get_value_or( T() ) + rhs.first.get_value_or( T() );
	if( lhs.second || rhs.second ) lhs.second = lhs.second.get_value_or( T() ) + rhs.second.get_value_or( T() );
	return lhs;
}


/** Implements x * log_y and handles special case where x=0 and log_y=+/-inf by returning 0. */
inline
double
safe_x_log_y( double x, double log_y ) {
	if( 0. == x ) {
		return 0.;
	} else {
		return x * log_y;
	}
}



/// Check the value is between 0 and 1.
inline
bool
is_between_0_and_1( double x ) {
	return 0. <= x && x <= 1.;
}



/// Find the index of the first unknown value or -1 if none
template< typename Seq >
inline
int
find_first_unknown( const Seq & seq ) {
	typedef typename seqan::Value< Seq >::Type value_t;
	int i = 0;
	BOOST_FOREACH( value_t b, seq ) {
		if( seqan::unknownValue< value_t >() == b ) {
			return i;
		}
		++i;
	}
	return -1;
}



/**
 * Functor that complements seqan::Dna5.
 */
struct ComplementDna5 : std::unary_function< seqan::Dna5, seqan::Dna5 > {

	/** Returns the complement of its argument. */
	seqan::Dna5 operator()( seqan::Dna5 b ) const {
		seqan::Dna5 result;
		switch( seqan::ordValue( b ) ) {
		case 0: result.value = 3; break;
		case 1: result.value = 2; break;
		case 2: result.value = 1; break;
		case 3: result.value = 0; break;
		case 4: result.value = 4; break;
		default: throw std::invalid_argument( "Unknown Dna5 base." );
		}
		return result;
	}
};


/**
 * Transforms an iterator into a reverse complement iterator.
 */
template< typename It >
boost::transform_iterator< ComplementDna5, boost::reverse_iterator< It > >
make_rev_comp_iterator( It it ) {
	return
		boost::make_transform_iterator(
			boost::make_reverse_iterator( it ),
			ComplementDna5()
		);
}


/**
 * Transforms a boost.range into its reverse complement.
 */
template< typename Range >
boost::iterator_range<
	boost::transform_iterator<
		ComplementDna5,
		boost::reverse_iterator<
			typename boost::range_iterator< const Range >::type
		>
	>
>
make_rev_comp_range( const Range & range ) {
	return
		boost::make_iterator_range(
			make_rev_comp_iterator( boost::end  ( range ) ),
			make_rev_comp_iterator( boost::begin( range ) )
		);
}





/// Default specification for STEME.
struct default_spec;



} // namespace steme



#endif /* STEME_JR_13AUG2011_DEFS_H_ */
