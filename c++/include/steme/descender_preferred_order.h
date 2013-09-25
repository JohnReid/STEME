/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the type that descends a suffix tree in a preferred order.
 *
 */

#ifndef STEME_JR_6SEP2011_DESCENDER_PREFERRED_ORDER_H_
#define STEME_JR_6SEP2011_DESCENDER_PREFERRED_ORDER_H_

#include <steme/defs.h>
#include <steme/descender.h>

#include <boost/tuple/tuple_comparison.hpp>

namespace steme {



/**
 * Knows how to descend a suffix tree choosing the preferred base at each level.
 */
template< typename Derived, typename Spec >
struct descender_preferred_order
: tree_descender< Derived, Spec >
{
    STEME_TYPEDEF_SEQAN_TYPES( Spec )
    typedef tree_descender< Derived, Spec >   descender_t;       ///< Descender base type.

    typedef std::vector< alphabet_t >         dna_vec;           ///< Vector of DNA type.
    typedef std::vector< dna_vec >            dna_matrix;        ///< Matrix of DNA type.

    dna_matrix                                preferred_bases;  ///< The preferred bases at each column of the model.

    /** Constructor. */
    descender_preferred_order(
        data< Spec > & _data,
        model< Spec > & _model
    )
    : descender_t( _data, _model )
    {
        preferred_bases.resize( _model.W() );
        for( size_t w = 0; size_t( _model.W() ) != w; ++w ) {
            calculate_preferred_bases_column( _model.bs.pssm.log_probs.m[ w ], preferred_bases[ w ] );
        }
    }

    /**
     * Calculate the preferred bases for a given column.
     */
    template< typename Col >
    static
    inline
    void
    calculate_preferred_bases_column( const Col & col, dna_vec & preferred ) {
        typedef boost::tuple< double, alphabet_t > sortable_t;
        typedef std::vector< sortable_t > to_sort_t;

        // build a vector that we can sort so most likely base is first
        to_sort_t to_sort;
        for( int i = 0; int( boost::size( col ) ) != i; ++i ) {
            to_sort.push_back( boost::make_tuple( col[ i ], alphabet_t( i ) ) );
        }
        std::sort( to_sort.begin(), to_sort.end(), std::greater< sortable_t >() );

        // copy the bases into the preferred bases vector.
        preferred.clear();
        for( unsigned int i = 0; boost::size( col ) != i; ++i ) {
            preferred.push_back( to_sort[ i ].template get< 1 >() );
        }
    }


    /**
     * Descend the suffix tree in a particular order.
     */
    inline
    void
    descend(
        top_down_it it,
        typename descender_t::partial_evaluation_pair partial_eval_pair,
        typename descender_t::bg_partial_evaluation_t partial_bg_eval,
        size_t w
    ) {
        MS_DEBUG_STRING_VAR( w_mer, representative( it ) );

        /// Visit rest of tree below node in preferred order.
        BOOST_FOREACH( alphabet_t base, preferred_bases[ repLength( it ) ] ) {
            //cout << "Descending: " << setw( W ) << representative( it ) << ": depth=" << rep_length << "; base=" << base << "\n";
            top_down_it i = it;
            if( goDown( i, base ) ) {
                //cout << "Went down from \"" << representative( it ) << "\" -> \"" << representative( i ) << "\"\n";
                descender_t::visit_node( i, partial_eval_pair, partial_bg_eval, w );
            }
        }
    }
};



} // namespace steme

#endif /* STEME_JR_6SEP2011_DESCENDER_PREFERRED_ORDER_H_ */
