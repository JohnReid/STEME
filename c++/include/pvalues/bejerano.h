/**
 * Copyright John Reid 2011
 *
 * \file Code that implements Bejerano's convex algorithm for exact p-value calculation.
 */


#ifndef _PVALUES_JR_22AUG2011_BEJERANO_H_
#define _PVALUES_JR_22AUG2011_BEJERANO_H_

#include <pvalues/defs.h>

namespace pvalues {



/**
 * Implementation of Bejerano's convex algorithm for exact p-value calculation.
 *
 * Notation similar to paper:
 *   Efficient Exact p-value Computation
 *   and Applications to Biosequence Analysis
 */
template< typename Spec = spec::Default >
struct bejerano {

    typedef typename spec::Float< Spec >::type float_t; ///< Floating point type.
    typedef typename spec::LogArithmetic< Spec >::type log_arithmetic_t;
    typedef typename log_arithmetic_t::log_t log_t; ///< Type that stores logarithms.
    typedef typename log_arithmetic_t::log_vec log_vec; ///< Vector of logarithms.
    typedef bejerano< Spec > self_t;
    typedef std::vector< float_t > float_vec; ///< Vector of doubles.
    typedef unsigned count_t; ///< Type that holds counts.
    typedef std::vector< count_t > count_vec; ///< Vector of counts.

    log_arithmetic_t & _log_arith; ///< Used to do log arithmetic.
    int _n; ///< Total observations.
    float_t _llr; ///< Log-likelihood ratio.
    float_t _llr_adj; ///< Adjusted LLR = LLR - n log(n)
    float_vec _q; ///< Background (null) distribution.
    log_vec _log_q; ///< Logarithms of background distribution.
    float_vec _q_bar; ///< Sum of rest of background distribution.
    log_vec _log_q_bar; ///< Logarithms of sum of rest of background distribution.
    float_vec _q_bar_min; ///< Sum of rest of background distribution.
    log_vec _log_q_bar_min; ///< Logarithms of sum of rest of background distribution.
    float_vec _q_bar_term; ///< q_i / (q_i + q_bar)
    float_vec _q_bar_min_term; ///< q_i / (q_i + q_bar_min)

    /// Construct.
    template< typename QRange >
    bejerano( log_arithmetic_t & log_arith, count_t n, float_t llr, const QRange & q )
    : _log_arith( log_arith )
    , _n( n )
    , _llr( llr )
    , _llr_adj( llr + n * log_arith.log( n ) )
    , _q( boost::begin( q ), boost::end( q ) )
    {
        using namespace boost::lambda;

        // calculate log(q)
        typename log_arithmetic_t::log_t ( log_arithmetic_t::* log_fn_ptr )( typename log_arithmetic_t::real_t ) = &log_arithmetic_t::log;
        std::copy(
            boost::make_transform_iterator( boost::begin( q ), bind( log_fn_ptr, &_log_arith, _1 ) ),
            boost::make_transform_iterator( boost::end( q ), bind( log_fn_ptr, &_log_arith, _1 ) ),
            std::back_inserter( _log_q )
        );
        if( spec::LogOn< Spec >::value ) {
            std::cout << "q:         ";
            std::copy( _q.begin(), _q.end(), std::ostream_iterator< float_t >( std::cout, " " ) );
            std::cout << "\n";
        }

        // calculate q_bar
        for( typename float_vec::iterator j = _q.begin(); ; ) {
            const float_t q = *j;
            const float_t q_bar = std::accumulate( ++j, _q.end(), 0. );
            if( _q.end() == j ) {
                break;
            }
            _q_bar.push_back( q_bar );
            _log_q_bar.push_back( _log_arith.log( q_bar ) );
            _q_bar_term.push_back( q / ( q + q_bar ) );
        }
        if( spec::LogOn< Spec >::value ) {
            std::cout << "q_bar:     ";
            std::copy( _q_bar.begin(), _q_bar.end(), std::ostream_iterator< float_t >( std::cout, " " ) );
            std::cout << "\n";
        }

        // calculate q_bar_min
        for( typename float_vec::iterator j = _q.begin(); ; ) {
            const float_t q = *j;
            typename float_vec::iterator k = std::min_element( ++j, _q.end() );
            if( _q.end() == j ) {
                break;
            }
            const float_t q_bar_min = *k;
            _q_bar_min.push_back( q_bar_min );
            _log_q_bar_min.push_back( _log_arith.log( q_bar_min ) );
            _q_bar_min_term.push_back( q / ( q + q_bar_min ) );
        }
        if( spec::LogOn< Spec >::value ) {
            std::cout << "q_bar_min: ";
            std::copy( _q_bar_min.begin(), _q_bar_min.end(), std::ostream_iterator< float_t >( std::cout, " " ) );
            std::cout << "\n";
        }
    }


    /// Calculate log p-value
    float_t
    operator()() {
        return descend( 0, _n, 0., _log_arith.log_fact( _n ) );
    }


protected:
    /// Descend tree of partial assignments
    log_t
    descend( count_t i, count_t n_bar, log_t llr_so_far, log_t log_Q_so_far ) const {
        PVALUE_ASSERT( n_bar );

        // result : log(p-value)
        log_t log_p_val = _log_arith.log_zero();

        // threshold on LLR
        const float_t bound_threshold = _llr_adj - llr_so_far;

        count_t alpha, beta, gamma, delta;

        // check if all sub-nodes are above threshold
        llr_bounds< false > llr_min( this, bound_threshold, i, n_bar );
        const double argmin_for_min = llr_min.argmin();
        if( ! llr_min( argmin_for_min ) ) {

            // all sub-nodes are above so we should add them all
            alpha = beta = gamma = delta = n_bar + 1;

        } else {

            // do two binary searches over LLR_min
            const count_t argmin_for_min_ceil = std::ceil( argmin_for_min );
            alpha = do_binary_search( llr_min, count_t( 0 ), argmin_for_min_ceil );
            delta = do_binary_search( make_unary_negate< count_t >( llr_min ), argmin_for_min_ceil, n_bar + 1 );

            // check if all maxima are above threshold
            llr_bounds< true > llr_max( this, bound_threshold, i, n_bar );
            const float_t argmin_for_max = llr_max.argmin();
            if( ! llr_max( argmin_for_max ) ) {

                // yes - so they must all be descended
                beta = gamma = delta;

            } else {

                // do two binary searches over LLR_max to determine beta and gamma
                const count_t argmin_for_max_ceil = std::ceil( argmin_for_max );
                beta = do_binary_search( llr_max, alpha, std::max( alpha, argmin_for_max_ceil ) );
                gamma = do_binary_search( make_unary_negate< count_t >( llr_max ), std::min( delta, argmin_for_max_ceil ), delta );

#ifndef PVALUE_CHECKS_DISABLED
                // check some invariants
                // const float_t argmin_for_min = llr_min.argmin();
                // const float_t argmin_for_max = llr_max.argmin();
                const float_t min_bound_at_argmin_for_min = llr_min.bound( argmin_for_min );
                const float_t min_bound_at_argmin_for_max = llr_min.bound( argmin_for_max );
                const float_t max_bound_at_argmin_for_min = llr_max.bound( argmin_for_min );
                const float_t max_bound_at_argmin_for_max = llr_max.bound( argmin_for_max );

                // check the ceil function worked
                PVALUE_ASSERT( argmin_for_min <= argmin_for_min_ceil );
                PVALUE_ASSERT( argmin_for_min_ceil - argmin_for_min < 1. );
                PVALUE_ASSERT( argmin_for_max <= argmin_for_max_ceil );
                PVALUE_ASSERT( argmin_for_max_ceil - argmin_for_max < 1. );

                // check bounds are not all above threshold, in that case we could have added this entire node in one go (we can still run into this with i = 0)
                PVALUE_ASSERT( i == 0 || llr_min( argmin_for_min ) );

                // check the minimums are better minimums for their bounds than the minimum for the other bound
                PVALUE_ASSERT( min_bound_at_argmin_for_min <= min_bound_at_argmin_for_max );
                PVALUE_ASSERT( max_bound_at_argmin_for_max <= max_bound_at_argmin_for_min );

                // check the real minimums are better minimums than the ceiled versions
                PVALUE_ASSERT( min_bound_at_argmin_for_min <= llr_min.bound( argmin_for_min_ceil ) );
                PVALUE_ASSERT( max_bound_at_argmin_for_max <= llr_max.bound( argmin_for_max_ceil ) );

                // check the max bound is bigger than the min bound at various points
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, 0u ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, alpha ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, beta ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, gamma ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, delta ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, n_bar ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, argmin_for_min ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, argmin_for_min_ceil ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, argmin_for_max ) );
                PVALUE_ASSERT( check_bounds( llr_min, llr_max, argmin_for_max_ceil ) );

                // check we found true cross-over points if we found any
                const float_t gamma_max = llr_max.bound( gamma );
                PVALUE_ASSERT( 0 == alpha || argmin_for_min_ceil == alpha || llr_min.bound( alpha ) <= bound_threshold );
                PVALUE_ASSERT( 0 == beta  || argmin_for_max_ceil == beta  || llr_max.bound( beta  ) <= bound_threshold );
                PVALUE_ASSERT( argmin_for_max_ceil == gamma || n_bar + 1 == gamma || llr_max.bound( gamma ) >= bound_threshold );
                PVALUE_ASSERT( argmin_for_min_ceil == delta || n_bar + 1 == delta || llr_min.bound( delta ) >= bound_threshold );

                // check some things are in the correct order
                PVALUE_ASSERT( 0 <= alpha );
                PVALUE_ASSERT( alpha <= beta );
                PVALUE_ASSERT( beta <= gamma );
                PVALUE_ASSERT( gamma <= delta );
                PVALUE_ASSERT( delta <= n_bar + 1 );
#endif // PVALUE_CHECKS_DISABLED

            }
        }

        // add the entire sub-trees that we can...
        for( count_t n_i = 0; alpha != n_i; ++n_i ) {
            log_p_val = _log_arith.log_add(
                log_p_val,
                subtree_log_Q( i, n_i, n_bar, log_Q_so_far )
            );
            PVALUE_ASSERT( _log_arith.log_isfinite( log_p_val ) );
        }
        for( count_t n_i = delta; n_bar + 1 != n_i; ++n_i ) {
            log_p_val = _log_arith.log_add(
                log_p_val,
                subtree_log_Q( i, n_i, n_bar, log_Q_so_far )
            );
            PVALUE_ASSERT( _log_arith.log_isfinite( log_p_val ) );
        }

        // descend the parts of the partial assignment tree that we need to
        for( count_t n_i = alpha; beta != n_i; ++n_i ) {
            log_p_val = _log_arith.log_add(
                log_p_val,
                descend_n_i( i, n_i, n_bar, llr_so_far, log_Q_so_far )
            );
            PVALUE_ASSERT( _log_arith.log_isfinite( log_p_val ) );
        }
        for( count_t n_i = gamma; delta != n_i; ++n_i ) {
            log_p_val = _log_arith.log_add(
                log_p_val,
                descend_n_i( i, n_i, n_bar, llr_so_far, log_Q_so_far )
            );
            PVALUE_ASSERT( _log_arith.log_isfinite( log_p_val ) );
        }

        return log_p_val;
    }

    inline
    float_t
    subtree_log_Q( count_t i, count_t n_i, count_t n_bar, log_t log_Q_so_far ) const {
        PVALUE_ASSERT( n_bar >= n_i );
        const count_t remaining = n_bar - n_i;
        const float_t result =
            log_Q_so_far
            + n_i * _log_q[ i ]
            - _log_arith.log_fact( n_i )
            + remaining * _log_q_bar[ i ]
            - _log_arith.log_fact( remaining )
            ;
        PVALUE_ASSERT( _log_arith.log_isfinite( result ) );
        PVALUE_ASSERT( result <= _log_arith.log( count_t( 1 ) ) );
        if( spec::LogOn< Spec >::value ) {
            std::cout << "subtree_log_Q(): "
                "i=" << i << "; "
                "n_i=" << n_i << "; "
                "n_bar=" << n_bar << "; "
                "log_Q_so_far=" << log_Q_so_far << "; "
                "result=" << result << "; "
                << "\n"
                ;
        }
        return result;
    }

    inline
    float_t
    descend_n_i( count_t i, count_t n_i, count_t n_bar, log_t llr_so_far, log_t log_Q_so_far ) const {
        const float_t log_q_i = _log_q[ i ];
        return descend(
            i + 1,
            n_bar - n_i,
            llr_so_far + n_i * ( _log_arith.log( n_i ) - log_q_i ),
            log_Q_so_far + n_i * log_q_i - _log_arith.log_fact( n_i )
        );
    }

    /**
     * Functional to calculate LLR_min or LLR_max, template parameter switches which log_q_bar_star is used.
     */
    template< bool Max >
    struct llr_bounds
    {
    protected:
        const self_t * _problem;
        float_t _threshold;
        count_t _i;
        count_t _n_bar;

    public:
        llr_bounds(
            const self_t * problem,
            float_t threshold,
            count_t i,
            count_t n_bar
        )
        : _problem( problem )
        , _threshold( threshold )
        , _i( i )
        , _n_bar( n_bar )
        { }

        /// \return The number left to be assigned.
        inline
        count_t
        n_bar() const {
            return _n_bar;
        }

        /// \return The position where the bound is minimised.
        inline
        float_t
        argmin() const {
            return ( Max ? _problem->_q_bar_min_term[ _i ] : _problem->_q_bar_term[ _i ] ) * _n_bar;
        }

        /// The bound.
        template< typename T >
        inline
        float_t
        bound( T n_i ) const {
            PVALUE_ASSERT( n_i <= _n_bar + 1 );
            const T remaining = _n_bar - n_i;
            float_t result = 0.;
            if( n_i ) { // 0 . log(0) == 0, so don't add in that case
                result += n_i * ( _problem->_log_arith.log( n_i ) - _problem->_log_q[ _i ] );
            }
            if( remaining ) { // 0 . log(0) == 0, so don't add in that case
                result += remaining * ( _problem->_log_arith.log( remaining ) - q_bar_star() );
            }
            PVALUE_ASSERT( boost::math::isfinite( result ) );
            return result;
        }

        /// Is the bound less than the threshold?
        template< typename T >
        bool
        operator()( T n_i ) const {
            const float_t b = bound( n_i );
            return b < _threshold;
        }

    protected:
        inline
        float_t
        q_bar_star() const {
            return Max ? _problem->_log_q_bar_min[ _i ] : _problem->_log_q_bar[ _i ] ;
        }

    };

    bool
    check_bounds( const llr_bounds< false > & llr_min, const llr_bounds< true > & llr_max, count_t n_i ) const {
        PVALUE_ASSERT( llr_min.n_bar() == llr_max.n_bar() );
        return
            n_i < 0u
            || llr_min.n_bar() < n_i
            || llr_min.bound( n_i ) <= llr_max.bound( n_i )
            ;
    }
};


} // namespace pvalues


#endif // _PVALUES_JR_22AUG2011_BEJERANO_H_
