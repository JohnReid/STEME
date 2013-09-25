/**
 * Copyright John Reid 2011
 *
 * \file Implementation of Nagarajan's shifted version of Hirji's algorithm.
 * See "Computing the P-value of the information content from an alignment of multiple
 * sequences" by Nagarajan, Jones and Keich
 */


#ifndef _PVALUES_JR_22AUG2011_HIRJI_SHIFTED_H_
#define _PVALUES_JR_22AUG2011_HIRJI_SHIFTED_H_

#include <pvalues/defs.h>

#define EULERS_NUMBER 2.71828182845904523536028747135266249775724709369995
#define ONE_OVER_EULERS_NUMBER ( 1. / EULERS_NUMBER )

namespace pvalues {




/**
 * Implementation of Nagarajan's shifted version of Hirji's algorithm.
 * See "Computing the P-value of the information content from an alignment of multiple
 * sequences" by Nagarajan, Jones and Keich
 */
template<
    typename Float = double
>
struct hirji_shifted {

    typedef hirji_shifted< Float > self_t; ///< This type.
    typedef std::vector< Float > float_vec; ///< Vector of floats.
    typedef std::vector< size_t > size_vec; ///< Vector of sizes.
    typedef std::vector< int > int_vec; ///< Vector of ints.
    typedef boost::multi_array< Float, 2 > array_t; ///< Array type.
    typedef boost::multi_array_ref< Float, 2 > array_ref_t; ///< Array reference type.
    typedef typename array_ref_t::const_reference array_ref_const_ref_t;
    typedef typename array_ref_const_ref_t::template const_array_view< 1 >::type slice_t;
    typedef typename array_t::reference array_row_ref; ///< Array row reference type.
    typedef boost::multi_array< int, 2 > int_array_t; ///< Integral array type.
    typedef typename int_array_t::reference int_array_row_ref; ///< Integral array row reference type.


    float_vec b; ///< Background distribution.
    const Float min_log_b; ///< Log of smallest background frequency.
    const size_t Q; ///< Size of lattice.
    const Float delta; ///< Delta that controls size of lattice
    const size_t n; ///< n we will calculate p.m.f. for
    array_t p_storage; ///< storage for p.m.f.s
    boost::scoped_ptr< array_ref_t > p; ///< p.m.f.s
    float_vec log_n; ///< log(n) for integral n.
    float_vec log_n_fact; ///< log(n!)
    int_array_t j; ///< j_a for each n'.
    array_t r; ///< r_a for each n'.
    size_vec argmin_j; ///< n1 at which minimum is achieved for each a, given delta and background frequencies.
    int_vec j_min; ///< Minimum possible j for each a, given delta and background frequencies.
    int_array_t p_first_j; /// Smallest j for which p_a_n > 0.
    int_array_t p_last_j; /// Largest j (+1) for which p_a_n > 0.


    /**
     * Holds cumulative distributions.
     */
    struct cumulative_distributions {
        typedef std::vector< float_vec > float_vec_vec; ///< Vector of vector of floats.

        float_t delta; ///< Lattice step size.
        float_vec adjustments; ///< Adjustments to LLR to fit lattice.
        float_vec_vec log_cmfs; ///< Log c.m.f.s

        cumulative_distributions( self_t & alg )
        : delta( alg.delta )
        , adjustments( alg.n )
        , log_cmfs( alg.n )
        {
            // build the c.m.f. for each n1
            BOOST_FOREACH( size_t n1, boost::counting_range( size_t( 1 ), alg.n + 1 ) ) {
                const double LLR_adjustment = alg.LLR_adjustment( n1 );
                adjustments[ n1 - 1 ] = LLR_adjustment - round( LLR_adjustment );
                const slice_t log_pmf = alg.get_log_pmf( n1 );
                double curr = std::log( 0. );
                log_cmfs[ n1 - 1 ].resize( boost::size( log_pmf ) );
                typename float_vec::reverse_iterator i = log_cmfs[ n1 - 1 ].rbegin();
                BOOST_FOREACH( float_t log_p, log_pmf | boost::adaptors::reversed ) {
                    curr = log_add( curr, log_p );
                    *i = curr;
                    ++i;
                }
            }
        }

        /// Round the LLR into a lattice index
        int
        round( float_t LLR ) const {
            return boost::math::round( LLR / delta );
        }

        /// Maximum n we have a c.m.f. for
        size_t
        max_n() const {
            return log_cmfs.size();
        }

        /// Get the log p-value.
        float_t
        get_log_p_value( float_t LLR, size_t n1 ) const {
            const size_t j = round( LLR );
            BOOST_ASSERT( n1 <= log_cmfs.size() );
            BOOST_ASSERT( j < log_cmfs[ n1 - 1 ].size() );
            if( n1 > log_cmfs.size() ) {
                throw std::invalid_argument(
                    PVALUE_MAKE_STRING(
                        "n1 is too big: n1=" << n1
                        << "; LLR=" << LLR
                        << "; delta=" << delta
                        << "; j=" << j
                        << "; log_cmfs.size()=" << log_cmfs.size()
                    )
                );
            }
            if( j >= log_cmfs[ n1 - 1 ].size() ) {
                throw std::invalid_argument(
                    PVALUE_MAKE_STRING(
                        "LLR is too big: LLR=" << LLR
                        << "; n1=" << n1
                        << "; delta=" << delta
                        << "; j=" << j
                        << "; log_cmfs[ n1 - 1 ].size()=" << log_cmfs[ n1 - 1 ].size()
                    )
                );
            }
            return log_cmfs[ n1 - 1 ][ j ];
        }
    };



    /// Construct.
    template< typename BgRange >
    hirji_shifted(
        const BgRange & bg,
        size_t n,
        size_t _Q = 0
    )
    : b( boost::begin( bg ), boost::end( bg ) )
    , min_log_b( boost::empty( bg ) ? 0. : std::log( *std::min_element( b.begin(), b.end() ) ) )
    , Q( _Q ? _Q : 2 * n )
    , delta( - ( n * min_log_b ) / ( Q - 1 ) )
    , n( n )
    , j( boost::extents[ boost::size( bg ) ][ n+1 ] )
    , r( boost::extents[ boost::size( bg ) ][ n+1 ] )
    , p_first_j( boost::extents[ boost::size( bg ) ][ n+1 ] )
    , p_last_j( boost::extents[ boost::size( bg ) ][ n+1 ] )
    {
        // needed to rebase our storage arrays
        typedef boost::multi_array_types::extent_range range;

        if( 0 == n ) {
            throw std::logic_error( "Cannot calculate p-values for n=0." );
        }
        if( 2 > boost::size( bg ) ) {
            throw std::logic_error( "Must have at least two background frequencies." );
        }

        // works quickest when background is sorted
        std::sort( b.begin(), b.end() );

        // calculate as much as possible up-front
        calculate_log_n();
        calculate_j();
        calculate_r();

        // storage
        const int minimal_j = get_minimal_j() - 1; // we can need this extra 1 later on due to rounding errors
        const int max_num_j = get_Q() - minimal_j;
        p_storage.resize( boost::extents[ n+1 ][ max_num_j ]);
        array_t tmp_p( boost::extents[ n+1 ][ max_num_j ]);
        const size_t num_elements = p_storage.num_elements();

        // calculate p for a = 0
        {
            const size_t a = 0;
            Float * data = ( A() - a ) % 2 ? p_storage.data() : tmp_p.data();
            std::fill( data, data + num_elements, Float( 0. ) );
            array_ref_t this_p(
                data,
                boost::extents[ n+1 ][ range( minimal_j, max_num_j + minimal_j ) ]
            );
            std::fill( p_first_j.data(), p_first_j.data() + p_first_j.num_elements(), std::numeric_limits< int >::max() );
            std::fill( p_last_j.data(), p_last_j.data() + p_last_j.num_elements(), - std::numeric_limits< int >::max() );
            int_array_row_ref j_a = j[ a ];
            for( size_t n1 = 0; n + 1 != n1; ++n1 ) {
                const int j_a_n1 = j_a[ n1 ];
                this_p[ n1 ][ j_a_n1 ] += r[ a ][ n1 ];
                p_first_j[ a ][ n1 ] = std::min( p_first_j[ a ][ n1 ], j_a_n1 );
                p_last_j[ a ][ n1 ] = std::max( p_last_j[ a ][ n1 ], j_a_n1 + 1 );
            }
        }

        // for all other a's
        for( size_t a = 1; a < A(); ++a ) {
            Float * data = ( A() - a ) % 2 ? p_storage.data() : tmp_p.data();
            Float * last_data = ( A() - a ) % 2 ? tmp_p.data() : p_storage.data();
            std::fill( data, data + num_elements, Float( 0. ) );
            array_ref_t this_p(
                data,
                boost::extents[ n+1 ][ range( minimal_j, max_num_j + minimal_j ) ]
            );
            array_ref_t last_p(
                last_data,
                boost::extents[ n+1 ][ range( minimal_j, max_num_j + minimal_j ) ]
            );
            int_array_row_ref j_a = j[ a ];
            array_row_ref r_a = r[ a ];
            for( size_t n1 = 0; n + 1 != n1; ++n1 ) {
                int smallest_j = std::numeric_limits< int >::max();
                int largest_j = -std::numeric_limits< int >::max();
                for( size_t n2 = 0; n1 + 1 != n2; ++n2 ) {
                    const size_t ndiff = n1 - n2;
                    const int j2 = j_a[ n2 ];
                    const double r2 = r_a[ n2 ];
                    const int first_j = p_first_j[ a - 1 ][ ndiff ] + j2;
                    const int last_j  = p_last_j[ a - 1 ][ ndiff ] + j2;
                    if( smallest_j > first_j ) { smallest_j = first_j; };
                    if( largest_j < last_j ) { largest_j = last_j; };
                    for( int j = first_j; last_j != j; ++j ) {
                        const int j_diff = j - j2;
                        BOOST_ASSERT( p_first_j[ a - 1 ][ ndiff ] <= j_diff );
                        BOOST_ASSERT( p_last_j[ a - 1 ][ ndiff ] > j_diff );
                        const Float lp = last_p[ ndiff ][ j_diff ];
                        if( lp ) {
                            this_p[ n1 ][ j ] += r2 * lp;
                            BOOST_ASSERT( boost::math::isfinite( this_p[ n1 ][ j ] ) );
                        }
                    }
                }
                p_first_j[ a ][ n1 ] = smallest_j; // remember smallest and largest for next time through loop
                p_last_j[ a ][ n1 ] = largest_j;
            }
        }

        // for final a
        {
            // get reference to last p
            const size_t a = A() - 1;
            Float * data = ( A() - a ) % 2 ? p_storage.data() : tmp_p.data();
            BOOST_ASSERT( p_storage.data() == data );
            p.reset(
                new array_ref_t(
                    data,
                    boost::extents[ n+1 ][ range( minimal_j, max_num_j + minimal_j ) ]
                )
            );

            for( size_t n1 = 0; n + 1 != n1; ++n1 ) {
                for( int j = minimal_j; max_num_j + minimal_j != j; ++j ) {
                    // scale p and make into log-value
                    ( *p )[ n1 ][ j ] = log_fact( n1 ) + std::log( ( *p )[ n1 ][ j ] ) - delta * j - n1 * ( log( n ) - 1. );

                    // must be log(0) or finite
                    BOOST_ASSERT( boost::math::isfinite( ( *p )[ n1 ][ j ] ) || std::log( 0. ) == ( *p )[ n1 ][ j ] );
                }
            }
        }
    }

    /// Smallest j for given n1
    inline
    int
    j_begin( size_t n1 = 0 ) const {
        if( ! n1 ) {
            n1 = n;
        }
        return p_first_j[ A() - 1 ][ n1 ];
    }

    /// One past maximum j for given n1
    inline
    int
    j_end( size_t n1 = 0 ) const {
        if( ! n1 ) {
            n1 = n;
        }
        return p_last_j[ A() - 1 ][ n1 ];
    }

    /// Value for given n1 and LLR
    inline
    Float
    get_p( Float LLR, size_t n1 = 0 ) const {
        return ( *p )[ n1 ][ get_j_for_LLR( LLR, n1 ) ];
    }

    /// j for given n1 and LLR
    inline
    int
    get_j_for_LLR( Float LLR, size_t n1 = 0 ) const {
        if( n1 ) {
            LLR += LLR_adjustment( n1 );
        }
        return round( LLR );
    }

    /// LLR for given n1 and j
    inline
    double
    get_LLR_for_j( int j, size_t n1 = 0 ) const {
        return j * delta - LLR_adjustment( n1 );
    }

    /// Value for given n1 and LLR
    inline
    Float
    get_p( int j, size_t n1 = 0 ) const {
        if( ! n1 ) {
            n1 = n;
        }
        return ( *p )[ n1 ][ j ];
    }

    /// Range for given n1
    inline
    slice_t
    get_log_pmf( size_t n1 = 0 ) {
        if( ! n1 ) {
            n1 = n;
        }

        // we can always achieve our highest j so just pick it out from pre-calculated bounds
        const int end_j = j_end( n1 );
        const int start_j = end_j + round( n1 * min_log_b ) - 1;

        // should have probability mass on last j
        BOOST_ASSERT( ( *p )[ n1 ][ end_j - 1 ] > log( 0 ) );

        // make sure we didn't go before first element
        BOOST_ASSERT( start_j >= p->index_bases()[ 1 ] );

        // due to rounding errors we may have probability mass in the lattice below the start so we move it
        // to the first position in the p.m.f.
        if( p_first_j[ A() - 1 ][ n1 ] < start_j ) {
            for( int j = start_j - 1; p_first_j[ A() - 1 ][ n1 ] <= j; --j ) {
                ( *p )[ n1 ][ start_j ] = log_add( ( *p )[ n1 ][ start_j ], ( *p )[ n1 ][ j ] );
                ( *p )[ n1 ][ j ] = std::log( 0. );
            }
            p_first_j[ A() - 1 ][ n1 ] = start_j;
        }

        typedef boost::multi_array_types::index_range range;
        return ( *p )[ n1 ][ boost::indices[ range( start_j, end_j ) ] ];
    }

    /// Maximum KL-divergence
    inline
    Float
    get_max_LLR( size_t n1 = 0 ) const {
        if( ! n1 ) {
            n1 = n;
        }
        return - ( n1 * min_log_b );
    }

    /// Size of lattice
    inline
    int
    get_Q() const {
        return Q;
    }

    /// Smallest j we can encounter in algorithm.
    inline
    int
    get_minimal_j() const {
        return std::min(
            std::accumulate( j_min.begin(), j_min.end(), 0 ),
            round( - ( n / EULERS_NUMBER ) )
        );
    }

    /// Size of alphabet.
    inline
    size_t
    A() const {
        return b.size();
    }



protected:

    /// Do rounding
    inline
    int
    round( Float LLR ) const
    {
        return boost::math::round( LLR / delta );
    }

    /// LLR adjustment for n1 != n
    inline
    Float
    LLR_adjustment( size_t n1 ) const {
        BOOST_ASSERT( n1 <= n );
        return n1 * ( log( n1 ) - log( n ) );
    }

    /// Log(n) for integral n
    inline
    Float
    log( size_t n1 ) const {
        return log_n[ n1 ];
    }

    /// Log(n!) for integral n
    inline
    Float
    log_fact( size_t n1 ) const {
        return log_n_fact[ n1 ];
    }

    /// Calculate log(n) up-front
    void
    calculate_log_n() {
        log_n.clear();
        log_n.reserve( n + 1 );
        log_n_fact.clear();
        log_n_fact.reserve( n + 1 );
        log_n.push_back( - std::numeric_limits< Float >::max() );
        log_n_fact.push_back( Float( 0. ) );
        for( size_t n1 = 0; n != n1; ++n1 ) {
            const Float lgn = std::log( n1 + 1 );
            BOOST_ASSERT( boost::math::isfinite( lgn ) );
            log_n.push_back( lgn );
            log_n_fact.push_back( lgn + *log_n_fact.rbegin() );
            BOOST_ASSERT( boost::math::isfinite( *log_n_fact.rbegin() ) );
        }
        BOOST_ASSERT( n + 1 == log_n.size() );
        BOOST_ASSERT( n + 1 == log_n_fact.size() );
    }

    /// Calculate j up-front
    void
    calculate_j() {

        // calculate all the j
        j_min.clear();
        for( size_t a = 0; A() != a; ++a ) {
            int_array_row_ref j_a = j[ a ];
            const double tmp = - log( n ) - std::log( b[ a ] );

            // work out the argmin and minimum
            const double argmin_j_a = n * b[ a ] * ONE_OVER_EULERS_NUMBER; // argmin when n is real
            const size_t argmin_j_a_floor = std::floor( argmin_j_a );
            const size_t argmin_j_a_ceil = std::ceil( argmin_j_a );
            const double min_j_a_n_floor = argmin_j_a_floor * ( log( argmin_j_a_floor ) + tmp );
            const double min_j_a_n_ceil = argmin_j_a_ceil * ( log( argmin_j_a_ceil ) + tmp );
            argmin_j.push_back( min_j_a_n_floor < min_j_a_n_ceil ? argmin_j_a_floor : argmin_j_a_ceil );
            const int min_j_a_n = round( ( min_j_a_n_floor < min_j_a_n_ceil ? min_j_a_n_floor : min_j_a_n_ceil ) );
            j_min.push_back( min_j_a_n );

            // calculate all the j's
            for( size_t n1 = 0; n + 1 != n1; ++n1 ) {
                const double unrounded_j = n1 * ( tmp + log( n1 ) );
                BOOST_ASSERT( boost::math::isfinite( unrounded_j ) );
                const int rounded_j = round( unrounded_j );
                j_a[ n1 ] = rounded_j;
                BOOST_ASSERT( min_j_a_n <= rounded_j );
                BOOST_ASSERT( argmin_j.back() != n1 || min_j_a_n == rounded_j );
            }
        }
    }

    /// \return Minimal j in range [0, n1]
    inline
    int
    get_min_j( size_t a, size_t n1 ) {
        // is n1 before the argmin?
        if( n1 < argmin_j[ a ] ) {
            return j[ a ][ n1 ]; // yes
        } else {
            return j_min[ a ];
        }
    }

    /// \return Maximal j in range [0, n1]
    inline
    int
    get_max_j( size_t a, size_t n1 ) {
        const int last = j[ a ][ n1 ]; // last j in range
        return last < 0 ? 0 : last; // we always have 0 in range, so is last bigger or not?
    }

    /// Calculate r up-front
    void
    calculate_r() {
        for( size_t a = 0; A() != a; ++a ) {
            int_array_row_ref j_a = j[ a ];
            array_row_ref r_a = r[ a ];
            const double tmp = log( n ) - 1. + std::log( b[ a ] );
            for( size_t n1 = 0; n + 1 != n1; ++n1 ) {
                r_a[ n1 ] = std::exp( n1 * tmp + delta * j_a[ n1 ] - log_fact( n1 ) );
                BOOST_ASSERT( boost::math::isfinite( r_a[ n1 ] ) );
            }
        }
    }

};



} // namespace pvalues

#endif // _PVALUES_JR_22AUG2011_HIRJI_SHIFTED_H_
