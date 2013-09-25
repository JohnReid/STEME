/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the type that finds instances of a motif.
 *
 */

#ifndef STEME_JR_6SEP2011_FIND_INSTANCES_H_
#define STEME_JR_6SEP2011_FIND_INSTANCES_H_

#include <steme/descender.h>
#include <steme/descender_preferred_order.h>
#include <steme/w_mer_evaluation.h>


namespace steme {


/**
 * Finds the instances of the model in the data which have Z above a given threshold.
 */
template< typename Spec = default_spec >
struct find_instances
: descender_preferred_order< find_instances< Spec >, Spec >
{

    STEME_TYPEDEF_SEQAN_TYPES( Spec )
    typedef find_instances< Spec >                                        this_t;            ///< This type.
    typedef boost::shared_ptr< this_t >                                   ptr;               ///< Pointer.
    typedef tree_descender< this_t, Spec >                                base_t;            ///< Base type.
    typedef w_mer_evaluation< Spec >                                      eval_t;            ///< Evaluation type.
    typedef typename eval_t::vec                                         instances_vec;     ///< Vector of instances.
    typedef typename eval_t::vec_ptr                                     instances_vec_ptr; ///< Pointer to instances.

    double                              Z_threshold;      ///< Z threshold at which to ignore instances.
    instances_vec_ptr                    instances;        ///< Instances.
    double                              llr_threshold;    ///< The log likelihood ratio threshold.


    /** Constructor. */
    find_instances(
        data< Spec > &         _data,
        model< Spec > &        _model,
        double                 Z_threshold
    )
    : descender_preferred_order< this_t, Spec >( _data, _model )
    , Z_threshold( Z_threshold )
    , instances( new instances_vec )
    , llr_threshold( std::log( Z_threshold / ( 1. - Z_threshold ) ) )
    {
        if( Z_threshold < 0 ) {
            throw std::invalid_argument( "Z threshold must be >= 0." );
        } else if( 1. < Z_threshold ) {
            throw std::invalid_argument( "Z threshold must be <= 1." );
        }
    }


    /**
     * Called by tree_descender base class as part of CRTP.
     */
    void
    on_descent_begin() {
        base_t::stats.reset();
        instances->clear();
    }


    /**
     * \return The threshold on log likelihood ratio
     */
    double
    get_llr_threshold() {
        return llr_threshold;
    }


    /// Handle a calculated Z for one strand of an occurrence of a W-mer
    template< bool RevComp >
    inline
    void
    handle_W_mer_occurrence_one_strand(
        global_pos_t global_pos,
        double Z_strand,
        double g_n,
        double lambda,
        double log_lambda,
        double log_p_X_n_theta_BS,
        double log_p_background
    ) {
        if( Z_strand >= Z_threshold ) {
            instances->push_back( w_mer_evaluation< Spec >( Z_strand, global_pos, RevComp ) );
        }
    }


    /// Handle all the W-mer occurrences of one node of the suffix tree
    inline
    void
    handle_W_mer_all_occurrences(
        top_down_it it,
        size_t num_occurrences,
        double Z_positive_total,
        double Z_negative_total
    ) {
    }
};






} // namespace steme

#endif /* STEME_JR_6SEP2011_FIND_INSTANCES_H_ */
