/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the model type for STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_MODEL_H_
#define STEME_JR_13AUG2011_MODEL_H_

#include <steme/data.h>
#include <steme/binding_model.h>
#include <steme/background_model.h>
#include <steme/w_mer_evaluation.h>

/**
 * Artificial limit on lambda. On real data sets this should never be used. It only comes into play
 * when the number of W-mers is less than or equal to the number of sites in the model.
 */
#define MAX_LAMBDA .99


namespace steme {




struct cached_lambda {

protected: // make sure when we update lambda, log_lambda_ratio is also updated
	double                     lambda;                         ///< Probability of binding site.
	double                     log_lambda;                     ///< log(lambda)
	double                     log_1_minus_lambda;             ///< log(1-lambda)
	double                     log_lambda_ratio;               ///< log_lambda - log_1_minus_lambda

	/// Update the logarithms of lambda and its inverse.
	void
	update_lambda_logs() {
		log_lambda = std::log( lambda );
		log_1_minus_lambda = std::log( 1. - lambda );
		log_lambda_ratio = log_lambda - log_1_minus_lambda;
	}

public:

	/// Constructor.
	cached_lambda( double lambda )
	: lambda( lambda )
	{
		update_lambda_logs();
	}

	/// \return lambda, the prior probability of a W-mer being a binding site.
	inline
	double
	get_lambda() const {
		return lambda;
	}

	/// \return log(lambda) - log(1 - lambda).
	inline
	double
	get_log_lambda_ratio() const {
		return log_lambda_ratio;
	}

	/// \return log(1 - lambda).
	inline
	double
	get_log_1_minus_lambda() const {
		return log_1_minus_lambda;
	}

	/// \return log(lambda).
	inline
	double
	get_log_lambda() const {
		return log_lambda;
	}

	/// set the prior lambda
	void
	set_lambda( double new_lambda ) {
		if( new_lambda > MAX_LAMBDA ) {
			new_lambda = MAX_LAMBDA;
		}
		lambda = new_lambda;
		update_lambda_logs();
	}
};





/// A model's prior.
template< typename Spec >
struct prior {
	cached_lambda                  lambda;             ///< Cached values of lambda
	typename psp< Spec >::ptr      _psp;               ///< Position-specific priors
	double                         num_sites;          ///< Prior estimate of number of sites.

	/// Constructor.
	prior( data< Spec > & _data, double lambda, double prior_num_sites, size_t W, bool dense_weights = true )
	: lambda( lambda )
	, _psp( _data.prior_for_width( W ) )
	, num_sites( prior_num_sites )
	{
	}

	/**
	 * Calculate the logarithm of the prior using the given g.
	 *
	 * The prior is lambda(1-g) / [1 - lambda(1-g)]
	 *
	 * Deals with the case that g==0 specially
	 */
	inline
	double
	calculate_log_prior( double g ) const {
        BOOST_ASSERT( 0 <= g );
        BOOST_ASSERT( g <= 1 );
	    if( g > 0. ) {
            const double prior_term = lambda.get_lambda() * ( 1. - g );
            return std::log( prior_term / ( 1. - prior_term ) );
	    } else {
	        return lambda.get_log_lambda_ratio();
	    }
	}

    /**
     * Calculate the the prior using the given g.
     */
    inline
    double
    calculate_prior( double g ) const {
        BOOST_ASSERT( 0 <= g );
        BOOST_ASSERT( g <= 1 );
        const double term = lambda.get_lambda() * ( 1. - g );
        return term / ( 1. - term );
    }
};



/**
 * Aggregates the components of the model : binding site, background, prior.
 */
template< typename Spec = default_spec >
struct model
: boost::noncopyable
{

	typedef boost::shared_ptr< model >                     ptr;                          ///< A pointer to a model.
    typedef data< Spec >                                  data_t;                       ///< The data type.
	typedef typename binding_model_meta< Spec >::type    bs_model_t;                   ///< The binding site model.
	typedef typename background_model_meta< Spec >::type bg_model_t;                   ///< The background model.
	STEME_TYPEDEF_SEQAN_TYPES( Spec )


    data_t &                     _data;              ///< The data.
	bs_model_t                   bs;                 ///< The binding site model.
	bg_model_t &                 bg;                 ///< The background model.
	prior< Spec >                _prior;             ///< Prior on binding sites.

protected:
	/// Protected copy constructor
	model( const model & rhs )
	: _data( rhs._data )
	, bs( rhs.bs )
	, bg( rhs.bg )
	, _prior( rhs._prior )
	{
	    // std::cout << "Copying model.\n";
	}

public:
	/** Constructor. If lambda out of range (0, 1) then assumes one site per sequence. */
	model(
		data_t & _data,
		bs_model_t bs,
		bg_model_t & bg,
		double lambda = 0.,
		double prior_num_sites = 0.
	)
		: _data( _data )
        , bs( bs )
		, bg( bg )
		, _prior( _data, lambda, prior_num_sites, W(), false )
	{
		if( lambda <= 0. || 1. <= lambda ) {
		    set_lambda_for_sites( _data.num_sequences() );
		}
	}

	/// Copy the model.
	ptr
	copy() {
		return ptr( new model( *this ) );
	}


	/** The width of the motif we are searching for. */
	size_t
	W() const {
		return bs.pssm.log_probs.W();
	}


	/** Calculate Z for the W-mer at the given position. Can over-ride the position-specific prior with optional argument. */
	double
	calculate_Z(
        global_pos_t global_pos,
        bool rev_comp,
        boost::optional< double > g = boost::optional< double >()
    ) {
	    // The position
        const int seq = seqan::getSeqNo( global_pos, _data.get_string_set_limits() );
        const int offset = seqan::getSeqOffset( global_pos, _data.get_string_set_limits() );

        // Likelihood of binding under model
	    const double log_p_binding =
            rev_comp
                ? w_mer_binding_log_likelihood< true  >( bs, _data.get_W_mer( W(), global_pos ) )
                : w_mer_binding_log_likelihood< false >( bs, _data.get_W_mer( W(), global_pos ) )
                ;

	    // Likelihood under background
	    const double log_p_background = bg.wmer_log_likelihood( seq, offset );

	    // Log Bayes factor
	    const double log_bayes_factor = log_p_binding - log_p_background;

	    // Prior
	    const double g_n = g ? *g : _prior._psp->get_g( seq, offset );
	    const double log_prior = _prior.calculate_log_prior( g_n );

	    // Calculate Z
	    return lr_to_prob( std::exp( log_prior + log_bayes_factor ) );
	}


	/**
	 * Get the W-mer that starts at the given index.
	 */
	typename seqan::Infix< string_t >::Type
	get_W_mer( int pos ) {
		return _data.get_W_mer( W(), pos );
	}


    /** Number of W-mers in index. */
    size_t
    num_W_mers() {
        return _data.num_W_mers( W() );
    }


    /** Number of known W-mers in index. */
    size_t
    num_occurrences() {
        return _data.get_occurrence_count( W() );
    }


	/** Set lambda for a given number of sites. */
	void
	set_lambda_for_sites( double num_sites ) {
        _prior.lambda.set_lambda( std::min( 1., num_sites / ( 2. * num_occurrences() ) ) );
	}

    /**
     * Updates the model with the best so many found sites. Also updates lambda.
     *
     * @return: How many sites used to update the model, this can be less than
     * requested if not enough were given.
     */
	template< typename Instances >
    size_t
    update_with_W_mers( const Instances & W_mers, size_t num_sites, bool use_pseudo_counts = false ) {

	    typename bs_model_t::new_model_t new_bs_model = bs.create_new_model( use_pseudo_counts );
        size_t sites_used = 0;
        typedef typename boost::range_value< Instances >::type w_mer_eval_t;
        BOOST_FOREACH( const w_mer_eval_t & W_mer, W_mers ) {
            if( sites_used == num_sites ) {
                break;
            }

            MS_DEBUG_STRING_RANGE( dbg_seq, get_W_mer( W_mer.global_pos ) );
            MS_DEBUG_STRING_RANGE( dbg_seq_rev_comp, make_rev_comp_range( get_W_mer( W_mer.global_pos ) ) );

            if( W_mer.rev_comp ) { // handle reverse complement case
                bs_model_t::update_new_model(
                    new_bs_model,
                    make_rev_comp_range( get_W_mer( W_mer.global_pos ) ), // W-mer
                    1. // weight
                );
            } else {
                bs_model_t::update_new_model(
                    new_bs_model,
                    get_W_mer( W_mer.global_pos ), // W-mer
                    1. // weight
                );
            }

            ++sites_used;
        }
        set_lambda_for_sites( sites_used );
        bs.update_from_new( new_bs_model ); //update the model using the newly learnt model

        return sites_used;
    }
};



} // namespace steme

#endif // STEME_JR_13AUG2011_MODEL_H_

