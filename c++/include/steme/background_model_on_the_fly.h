/** Copyright John Reid 2011, 2012
 *
 * \file
 * \brief Defines the on the fly background model type for STEME algorithm.
 *
 */

#ifndef STEME_JR_11JUL2011_BACKGROUND_MODEL_ON_THE_FLY_H_
#define STEME_JR_11JUL2011_BACKGROUND_MODEL_ON_THE_FLY_H_

#include <steme/seqan_types.h>
#include <steme/background_model.h>


namespace steme {




/// Calculates the likelihood of a infix in context on the fly.
template< typename Spec = default_spec >
struct on_the_fly_evaluator {
    virtual ~on_the_fly_evaluator() { }

    typedef boost::shared_ptr< on_the_fly_evaluator< Spec > >  ptr;  ///< Pointer.
    STEME_TYPEDEF_SEQAN_TYPES( Spec )
    typedef typename seqan::Iterator< prefix_t >::Type it;

    virtual double evaluate( it context, it begin, it end ) = 0;
    virtual double min_ll( size_t w, size_t W ) = 0;
};




/**
 * A background model for a particular motif width, W. This model calculates the
 * log likelihoods on the fly. This may be slightly slower but minimises the storage requirements so that
 * it is suitable when the indexed text is large, e.g. a genome.
 */
template< typename Spec = default_spec >
struct on_the_fly_bg_model
{
    typedef boost::shared_ptr< on_the_fly_bg_model >  ptr;  ///< Pointer.
    typedef double partial_evaluation_t; ///< A partial evaluation
    typedef typename on_the_fly_evaluator< Spec >::ptr evaluator_t; ///< An evaluator
    STEME_TYPEDEF_SEQAN_TYPES( Spec )


    evaluator_t evaluator;                           ///< Markov model of the background.


    on_the_fly_bg_model( evaluator_t evaluator )
    : evaluator( evaluator )
    {
    }


    /// Return a null partial evaluation.
    static
    partial_evaluation_t
    empty_partial_evaluation() {
        return 0.;
    }


    /// Update the partial evaluation.
    inline
    partial_evaluation_t
    update_partial_eval( prefix_t pre, partial_evaluation_t partial_eval, size_t w ) {
        return partial_eval + evaluator->evaluate( seqan::begin( pre ), seqan::begin( pre ) + w, seqan::end( pre ) );
    }


    /// Return the minimum log likelihood of any W-mer that starts with the prefix represented by this iterator.
    inline
    double
    get_min_log_likelihood( top_down_it, partial_evaluation_t partial_eval, size_t w, size_t W ) const {
        return partial_eval + evaluator->min_ll( w, W );
    }

    inline
    double
    wmer_log_likelihood( size_t seq, size_t offset ) const {
        throw std::runtime_error( "wmer_log_likelihood(): not implemented.");
    }

    /// Do we have log likelihoods at the base pair level?
    inline
    bool
    have_base_pair_LLs() {
        return false; // no
    }
};







} //namespace steme

#endif //STEME_JR_11JUL2011_BACKGROUND_MODEL_ON_THE_FLY_H_
