/** Copyright John Reid 2011, 2012
 *
 * \file
 * \brief Defines the suffix tree descender type for STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_DESCENDER_H_
#define STEME_JR_13AUG2011_DESCENDER_H_

#include <steme/data.h>
#include <steme/model.h>
#include <boost/timer.hpp>

namespace steme {



/**
 * Used to keep track of how many nodes in the suffix tree we are discarding/evaluating.
 */
struct efficiency_stats : boost::addable< efficiency_stats > {

	/** Counts how many discarded/evaluated at a depth. */
	struct count
	: boost::addable< count >
	, boost::less_than_comparable< count >
	, boost::equality_comparable< count >
	{
		typedef std::vector< count >   vec;          ///< Vector of counts.

		size_t                         evaluated;    ///< # evaluated at this depth.
		size_t                         discarded;    ///< # evaluated at this depth.

		/// Constructor.
		count() : evaluated( 0 ), discarded( 0 ) { }

		/// Reset counts.
		void reset() { evaluated = discarded = 0; }

		/// Total
		size_t total() const { return evaluated + discarded; }

		/// Add count to this one.
		count & operator+=( const count & rhs ) {
			this->evaluated += rhs.evaluated;
			this->discarded += rhs.discarded;
			return *this;
		}

		/// Equality test
		bool operator==( const count & rhs ) const {
			return this->evaluated == rhs.evaluated && this->discarded == rhs.discarded;
		}

		/// Comparison test
		bool operator<( const count & rhs ) const {
			if( this->fraction_evaluated() < rhs.fraction_evaluated() ) return true;
			if( this->fraction_evaluated() > rhs.fraction_evaluated() ) return false;
			return this->evaluated < rhs.evaluated;
		}

		/// Fraction evaluated
		double fraction_evaluated() const {
			return double( evaluated ) / (evaluated + discarded);
		}
	};




	/** Counts how many nodes & occurrences evaluated/discarded at a depth. */
	struct node_counts
	: boost::addable< node_counts >
	, boost::less_than_comparable< node_counts >
	, boost::equality_comparable< node_counts >
	{
		typedef std::vector< node_counts >   vec;       ///< Vector of counts.

		count         node_count;                       ///< Node counts.
		count         occurrence_count;                 ///< Occurrence counts.

		/// Reset counts.
		void reset() {
			occurrence_count.reset();
			node_count.reset();
		}

		/// Add count to this one.
		node_counts & operator+=( const node_counts & rhs ) {
			this->node_count += rhs.node_count;
			this->occurrence_count += rhs.occurrence_count;
			return *this;
		}

		/// Equality test
		bool operator==( const node_counts & rhs ) const {
			return this->node_count == rhs.node_count && this->occurrence_count == rhs.occurrence_count;
		}

		/// Comparison test
		bool operator<( const node_counts & rhs ) const {
			if( this->node_count < rhs.node_count ) return true;
			if( this->node_count > rhs.node_count ) return false;
			return this->occurrence_count < rhs.occurrence_count;
		}
	};

	node_counts::vec         counts;             ///< # nodes evaluated/discarded at each depth.

	/** Constructor. */
	efficiency_stats( size_t max_w = 0 )
	: counts( max_w + 1 )
	{
	}

	/** Add another set of efficiency statistics to this one. */
	efficiency_stats & operator+=( const efficiency_stats & rhs ) {
		if( counts.size() != rhs.counts.size() ) {
			throw std::logic_error( "Cannot add efficiency statistics of differing lengths" );
		}
		for( size_t i = 0; counts.size() != i; ++i ) {
			counts[ i ] += rhs.counts[ i ];
		}
		return *this;
	}

	/** Resets the statistics. */
	void
	reset() {
		for( size_t i = 0; counts.size() != i; ++i ) {
			counts[ i ].reset();
		}
	}

	/** Called when discarding a node. */
	void
	discarded( size_t w, size_t num_occurrences ) {
		MYRRH_ASSERT( w < counts.size() );
		counts[ w ].node_count.discarded++;
		counts[ w ].occurrence_count.discarded += num_occurrences;
	}

	/** Called when we don't discard a node. */
	void
	evaluated( size_t w, size_t num_occurrences ) {
		MYRRH_ASSERT( w < counts.size() );
		counts[ w ].node_count.evaluated++;
		counts[ w ].occurrence_count.evaluated += num_occurrences;
	}
};



/**
 * Descend the suffix tree estimating bounds as we go. Used as base class in CRTP idiom. The derived
 * class implements a threshold function that allows the tree descender to ignore some of the
 * branches.
 */
template<
	typename Derived,
	typename Spec = default_spec
>
struct tree_descender
{

	STEME_TYPEDEF_SEQAN_TYPES( Spec )
	typedef data< Spec >                                                 data_t;                    ///< Data type.
	typedef model< Spec >                                                model_t;                   ///< Model type.
    typedef typename model_t::bs_model_t                                bs_model_t;                ///< Binding site model type.
    typedef typename model_t::bg_model_t                                bg_model_t;                ///< Background model type.
    typedef typename model_t::bs_model_t::partial_evaluation_t          bs_partial_evaluation_t;   ///< The type that stores partial evaluations by the binding site model.
    typedef typename model_t::bg_model_t::partial_evaluation_t          bg_partial_evaluation_t;   ///< The type that stores partial evaluations by the background model.
    typedef optional_pair< bs_partial_evaluation_t >                     partial_evaluation_pair;   ///< A pair of optional partial evaluations by the binding site model.

	data_t &               _data;                   ///< The data.
	model_t &              _model;                  ///< The model of our data.
	efficiency_stats       stats;                   ///< Statistics on how many nodes we dropped etc...
	size_t 	               W;                       ///< The width of the model.


	/** Constructor. */
	tree_descender(
		data_t &                 _data,
		model_t &                _model
	)
	: _data( _data )
	, _model( _model )
	, stats( _model.W() )
	, W( _model.W() )
	{ }


	/**
	 * Alternate implementation to descend tree.
	 */
	void
	descend_tree() {
		static_cast< Derived * >( this )->on_descent_begin();

#ifdef STEME_LOG_TIMING_INFO
        boost::timer t;
#endif // STEME_LOG_TIMING_INFO

		visit_node(
			top_down_it( _data.index ),
			partial_evaluation_pair( bs_model_t::empty_partial_evaluation(), bs_model_t::empty_partial_evaluation() ),
			bg_model_t::empty_partial_evaluation(),
			0
		);

#ifdef STEME_LOG_TIMING_INFO
        std::cout << "descend_tree() took " << t.elapsed() << " seconds.\n";
#endif // STEME_LOG_TIMING_INFO

		static_cast< Derived * >( this )->on_descent_end();
	}


	/**
	 * Alternate implementation to descend tree.
	 */
	void
	visit_node(
		top_down_it it,
		partial_evaluation_pair partial_eval_pair,
		bg_partial_evaluation_t partial_bg_eval,
		size_t previous_w
	) {
		using namespace seqan;

		// we must have an eval on at least one strand if we've got here.
		MYRRH_ASSERT( partial_eval_pair );

		MS_DEBUG_STRING_VAR( w_mer, representative( it ) );

		//ignore if N in parent edge
		if( _data.contains_n( W, it ) ) {
			return;
		}

        // get bounds on g
        const double min_g = _model._prior._psp->get_min_g( it );

        // the depth we're at
        const size_t w = std::min( repLength( it ), W );
        MYRRH_ASSERT( w > previous_w || 0 == previous_w );

        // the number of possible binding sites
        const size_t num_occurrences = countOccurrences( it );

        // if g = 1. then we cannot have a binding site down this part of the suffix tree
        if ( min_g >= 1. ) {

            stats.discarded( w, num_occurrences ); // update stats

        } else {

            // the prefix represented by this node.
            const prefix_t pre = prefix( representative( it ), w );

            // update our partial evaluation of the background model
            partial_bg_eval = _model.bg.update_partial_eval( pre, partial_bg_eval, previous_w );

            // get a bound on the evaluation under the background model
            const double bg_log_likelihood_bound = _model.bg.get_min_log_likelihood( it, partial_bg_eval, w, W );

            // the threshold on log likelihood ratio
            const double llr_threshold = static_cast< Derived * >( this )->get_llr_threshold();

            // log prior
            const double log_prior = _model._prior.calculate_log_prior( min_g );

            // get the bound on log p(X_n|\theta_BS)
            const double eval_threshold =
                llr_threshold
                - log_prior
                + bg_log_likelihood_bound
                ;

            // we have to evaluate the W-mer and its reverse complement.
            // start with the W-mer
            if( partial_eval_pair.first ) { // do we have an eval on the positive strand?
                partial_eval_pair.first =
                    _model.bs.template update_partial_eval_thresholded< false >(
                        pre,
                        *partial_eval_pair.first,
                        previous_w,
                        eval_threshold
                    );
            }

            // now evaluate the reverse complement
            if( partial_eval_pair.second ) { // do we have an eval on the negative strand?
                partial_eval_pair.second =
                    _model.bs.template update_partial_eval_thresholded< true >(
                        pre,
                        *partial_eval_pair.second,
                        previous_w,
                        eval_threshold
                    );
            }

            // do we have any evals left?
            if( ! partial_eval_pair ) { // do we have an eval on either strand?

                stats.discarded( w, num_occurrences ); // update stats

            } else {

                stats.evaluated( w, num_occurrences ); // update stats

                // did we get to the bottom?
                if( W == w ) {

                    // yes - so evaluate the W-mer
                    static_cast< Derived * >( this )->evaluate_W_mer(
                        it,
                        partial_eval_pair,
                        bg_log_likelihood_bound,
                        num_occurrences
                    );

                } else {

                    // can be over-ridden in derived class.
                    static_cast< Derived * >( this )->descend( it, partial_eval_pair, partial_bg_eval, w );

                }
            }
        }
	}


	/**
	 * Descend the suffix tree. Can be over-ridden in derived classes.
	 */
	inline
	void
	descend(
		top_down_it it,
		partial_evaluation_pair partial_eval_pair,
        bg_partial_evaluation_t partial_bg_eval,
		size_t w
	) {
		// continue descending the suffix tree, so go down and across
		if( goDown( it ) ) {
			visit_node( it, partial_eval_pair, partial_bg_eval, w );
			while( goRight( it ) ) {
				visit_node( it, partial_eval_pair, partial_bg_eval, w );
			}
		}
	}


	/**
	 * Default implementation to be called as part of CRTP.
	 */
	void
	on_descent_begin() {
	}


	/**
	 * Default implementation to be called as part of CRTP.
	 */
	void
	on_descent_end() {
	}


	/**
	 * Default implementation to be called as part of CRTP.
	 * Calls evaluate_W_mer_strand in the derived class for each
	 * strand.
	 *
	 * Can be over-ridden in base classes.
	 */
	GCC_DIAG_OFF(GCC_UNINIT_WARNING); // turn off spurious uninit warnings
	void
	evaluate_W_mer(
		top_down_it it,
		partial_evaluation_pair bs_eval,
		double log_p_background,
		size_t num_occurrences
	) {
		MS_DEBUG_STRING_VAR( w_mer, representative( it ) );
		MYRRH_ASSERT( bs_eval ); // must have at least one eval here.

		const double lambda = _model._prior.lambda.get_lambda();
		const double log_lambda = _model._prior.lambda.get_log_lambda();

		// calculate likelihoods of X_n given binding model and log Bayes factors
		double log_p_X_n_theta_BS_positive, lbf_positive;
		if( bs_eval.first ) {
			log_p_X_n_theta_BS_positive = _model.bs.eval_to_log( *bs_eval.first );
            if( ! _model.bg.have_base_pair_LLs() ) {
                lbf_positive = log_p_X_n_theta_BS_positive - log_p_background;
            }
		}
		double log_p_X_n_theta_BS_negative, lbf_negative;
		if( bs_eval.second ) {
			log_p_X_n_theta_BS_negative = _model.bs.eval_to_log( *bs_eval.second );
			if( ! _model.bg.have_base_pair_LLs() ) {
			    lbf_negative = log_p_X_n_theta_BS_negative - log_p_background;
			}
		}

		// total Zs over all occurrences
		double Z_positive_total = 0.;
		double Z_negative_total = 0.;

		// for each occurrence, store its Z so it can be normalised later.
		occurrences_t occs = seqan::getOccurrences( it );
		MYRRH_ASSERT( seqan::length( occs ) == num_occurrences );
		for( size_t i = 0; num_occurrences != i; ++i ) {
			// update Z and evaluation vectors

			// get this occurrence
			const typename seqan::Value< occurrences_t >::Type & occ = occs[ i ];

			// get the global position to know where to store the Z
			const global_pos_t global_pos = seqan::posGlobalize( occ, _data.get_string_set_limits() );

			// prior
			const double g_n = _model._prior._psp->get_g( occ.i1, occ.i2 );
			const double log_prior = _model._prior.calculate_log_prior( g_n );

			// Z for the positive strand.
			if( bs_eval.first ) {
			    // if we are using a different LL for each occurrence
	            if( _model.bg.have_base_pair_LLs() ) {
	                lbf_positive = log_p_X_n_theta_BS_positive - _model.bg.wmer_log_likelihood( occ.i1, occ.i2 );
	            }
				const double Z = lr_to_prob( std::exp( lbf_positive + log_prior ) );
				MS_CHECK_BOUNDS( Z, "Z_positive", 0., 1. );
				Z_positive_total += Z;
				static_cast< Derived * >( this )->template handle_W_mer_occurrence_one_strand< false >(
					global_pos,
					Z,
					g_n,
					lambda,
					log_lambda,
					log_p_X_n_theta_BS_positive,
					log_p_background
				);
			}

			// Z for the negative strand.
			if( bs_eval.second ) {
                // if we are using a different LL for each occurrence
                if( _model.bg.have_base_pair_LLs() ) {
                    lbf_negative = log_p_X_n_theta_BS_negative - _model.bg.wmer_log_likelihood( occ.i1, occ.i2 );
                }
                const double Z = lr_to_prob( std::exp( lbf_negative + log_prior ) );
				MS_CHECK_BOUNDS( Z, "Z_negative", 0., 1. );
				Z_negative_total += Z;
				static_cast< Derived * >( this )->template handle_W_mer_occurrence_one_strand< true >(
					global_pos,
					Z,
					g_n,
					lambda,
					log_lambda,
					log_p_X_n_theta_BS_negative,
					log_p_background
				);
			}
		}

		// handle all W-mer occurrences
		static_cast< Derived * >( this )->handle_W_mer_all_occurrences(
			it,
			num_occurrences,
			Z_positive_total,
			Z_negative_total
		);
	}
    GCC_DIAG_ON(GCC_UNINIT_WARNING);
};







} // namespace steme

#endif /* STEME_JR_13AUG2011_DESCENDER_H_ */
