/** Copyright John Reid 2011, 2012
 *
 * \file
 * \brief Defines the expectation-maximization descender type for STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_EM_H_
#define STEME_JR_13AUG2011_EM_H_

#include <steme/descender.h>


namespace steme {



/**
 * Descends the suffix tree performing EM algorithm.
 *
 * \todo Revise code that checks bound on approximation
 *
 * \todo Move from storing a separate evaluation for each strand to just one evaluation? This would need
 * major work.
 */
template< typename Spec = default_spec >
struct EM_descender
: tree_descender< EM_descender< Spec >, Spec >
{

	STEME_TYPEDEF_SEQAN_TYPES( Spec )
	typedef tree_descender< EM_descender< Spec >, Spec >           base_t;          ///< Base type.
	typedef boost::shared_ptr< EM_descender< Spec > >              ptr;             ///< Shared pointer.
	typedef typename base_t::bs_model_t::new_model_t               new_model_t;     ///< The type of the new models we update in EM.
	typedef optional_pair< double >                                Z_pair_t;        ///< Represents optional evaluations (Zs) on the positive and negative strands.
	typedef boost::numeric::ublas::mapped_vector< Z_pair_t >       Z_sparse_vec;    ///< A sparse vector of evaluations.
	//typedef boost::numeric::ublas::compressed_vector< Z_pair_t >   Z_vec;         ///< A sparse vector of evaluations.
	//typedef boost::numeric::ublas::coordinate_vector< Z_pair_t >   Z_vec;         ///< A sparse vector of evaluations.
	typedef typename Z_sparse_vec::difference_type                 Z_sparse_diff_t; ///< Difference type for the Z vector.
	typedef typename Z_sparse_vec::size_type                       Z_sparse_size_t; ///< Size type for the Z vector.
	typedef typename Z_sparse_vec::const_iterator                  Z_sparse_it;     ///< Iterator into Z.
	typedef std::vector< Z_pair_t >                                Z_dense_vec;     ///< A dense vector of evaluations.
	typedef typename Z_dense_vec::difference_type                  Z_dense_diff_t;  ///< Difference type for the Z vector.
	typedef typename Z_dense_vec::size_type                        Z_dense_size_t;  ///< Size type for the Z vector.
	typedef typename Z_dense_vec::iterator                         Z_dense_it;      ///< Iterator into Z.
	typedef std::vector< int >                                     index_vec;       ///< Vector to hold indices into Z.


	double                     epsilon;                  ///< Allowed relative error in f_wb.
	double                     estimated_c_wb;           ///< The estimate of c_wb we used at the beginning of the iteration.
	double                     log_delta_factor;         ///< The amount by which we change log delta when we see the estimation of c_wb is too low.
	double                     log_delta_adjustment_max; ///< The largest adjustment to log delta we allow.
	double                     log_delta_adjustment;     ///< An adjustment we make to log delta for the next iteration after seeing the actual value of c_wb. This can change every iteration.
	double                     llr_threshold;            ///< Threshold on log likelihood ratio.
	double                     wnsites;                  ///< Weight on number of sites. Used when updating lambda.
	bool                       using_sparse_Z;           ///< Are we using a sparse Z?
	Z_sparse_vec               Z_sparse;                 ///< The evaluations (expected Zs) of each W-mer under the binding site model.
	Z_dense_vec                Z_dense;                  ///< The evaluations (expected Zs) of each W-mer under the binding site model.
	double                     large_Z_threshold;        ///< The threshold at which Z are large (might need normalising).
	index_vec                  large_Z_indices;          ///< The indices of the Z over a threshold.
	double                     expected_sites;           ///< The expected number of sites.
	double                     LL;                       ///< Expected log likelihood.
	double                     log_lambda_ratio;         ///< log(lambda) - log(1-lambda).
	std::vector< double >      lambda_ratios;            ///< Store the lambda ratios. These tell us whether our threshold was high enough.
	new_model_t                new_bs_model;             ///< Holds the new model we are updating.
	size_t                     num_Z_large;              ///< The number of Z that were large.
	size_t                     num_Z_non_zero;           ///< The number of Z that were non-zero.
	size_t                     num_Z_normalised;         ///< The number of windows that were normalised.


	/** Constructor. */
	EM_descender(
		typename base_t::data_t &       _data,
		typename base_t::model_t &      _model,
		double                          epsilon,
		double                          wnsites,
		double                          log_delta_factor = std::log( .9 ),
		double                          log_delta_adjustment_max = std::log( 5. )
	)
	: base_t( _data, _model )
	, epsilon( epsilon )
	, estimated_c_wb( 0. )
	, log_delta_factor( log_delta_factor )
	, log_delta_adjustment_max( log_delta_adjustment_max )
	, log_delta_adjustment( 0. )
	, llr_threshold( -std::numeric_limits< double >::max() )
	, wnsites( wnsites )
	, using_sparse_Z( false )
	, Z_sparse( _data.N, _data.N / 3 )
	, Z_dense( _data.N )
	, large_Z_threshold( 1. / base_t::W )
	, expected_sites( 0 )
	, LL( 0. )
    , log_lambda_ratio( base_t::_model._prior.lambda.get_log_lambda_ratio() )
    , num_Z_large( 0 )
    , num_Z_non_zero( 0 )
    , num_Z_normalised( 0 )
	{
	}


	/** Get the Z associated with index. */
	Z_pair_t
	get_Z( size_t n ) const {
		return using_sparse_Z ? Z_sparse[ n ] : Z_dense[ n ];
	}


	/** Do one iteration. */
	double
	do_iteration() {
	    base_t::descend_tree();
		return expected_sites;
	}

	/**
	 * Called by tree_descender base class as part of CRTP.
	 */
	void
	on_descent_begin() {
		// set the threshold
	    const double delta_threshold = epsilon * base_t::_model._prior.lambda.get_lambda();
		llr_threshold = std::log( delta_threshold ) - std::log( 1. - delta_threshold );
		estimated_c_wb = base_t::_model._prior.lambda.get_lambda() * base_t::_data.get_occurrence_count( base_t::_model.W() );

		//reset Z estimates
		if( using_sparse_Z ) {
			Z_sparse.clear();
		} else {
			Z_dense.clear();
			Z_dense.resize( base_t::_data.N, Z_pair_t() );
		}
		MYRRH_ASSERT(
			( using_sparse_Z && Z_sparse.size() == base_t::_data.N )
			|| ( ! using_sparse_Z && Z_dense.size() == base_t::_data.N )
		);

		// clear the large Z index vector
		large_Z_indices.clear();

		// reset counts
		num_Z_large = 0;
		num_Z_non_zero = 0;
		num_Z_normalised = 0;

		// create a new model to update with the weighted sites.
		new_bs_model = base_t::_model.bs.create_new_model( true );

		// initialise how many sites the model expects
		expected_sites = 0.0;

		// initialise the stuff for the LL
		log_lambda_ratio = base_t::_model._prior.lambda.get_log_lambda_ratio();
		LL = 0.; //reset log likelihood
	}


	/**
	 * Called by tree_descender base class as part of CRTP.
	 */
	void
	on_descent_end() {
		// reset Z vectors we are not using
		if( using_sparse_Z ) {
			Z_dense.clear();
			shrink_to_empty( Z_dense );
		} else {
			Z_sparse.clear();
		}
		MYRRH_ASSERT(
			( using_sparse_Z && Z_sparse.size() == base_t::_data.N )
			|| ( ! using_sparse_Z && Z_dense.size() == base_t::_data.N )
		);

		// normalise those neighbouring Z that sum to more than 1.
		normalise_Z();

		// make adjustment if we over-estimated c_wb to start with.
		// this will be used on the next iteration
		if( expected_sites < estimated_c_wb ) {
			//std::cout << "WARNING: expected_sites=" << expected_sites << " < " << estimated_c_wb << "=estimated_c_wb: Relative error bounds may not be satisfied!\n";
			log_delta_adjustment += log_delta_factor;
		} else {
			log_delta_adjustment -= log_delta_factor;
			log_delta_adjustment = std::min( log_delta_adjustment, log_delta_adjustment_max ); // don't let it go too high
		}

		// update the model with the new model.
		base_t::_model.bs.update_from_new( new_bs_model ); //update the model using the newly learnt model

		//adjust num_sites by weight
		double adjusted_num_sites = expected_sites * (1 - wnsites) + base_t::_model._prior.num_sites * wnsites;

		// update the model with the adjusted lambda
		base_t::_model.set_lambda_for_sites( adjusted_num_sites );
		lambda_ratios.push_back( std::exp( base_t::_model._prior.lambda.get_log_lambda_ratio() ) );

//			const size_t num_W_mers = 2 * base_t::_data.num_W_mers( base_t::W );
//			std::cout
//				<< "# W-mers = " << num_W_mers << "; "
//				<< "# non-zero = " << std::setw( 8 ) << num_Z_non_zero << "=" << std::setw( 3 ) << int( 100 * num_Z_non_zero / num_W_mers ) << "%; "
//				<< "# above threshold = " << std::setw( 8 ) << num_Z_large << "=" << std::setw( 3 ) << int( 100 * num_Z_large / num_Z_non_zero ) << "%; "
//				<< "# normalised = " << std::setw( 8 ) << num_Z_normalised << "=" << std::setw( 3 ) << int( 100 * num_Z_normalised / num_Z_non_zero ) << "%"
//				<< "\n"
//				;
	}


	/**
	 * \return The threshold on log likelihood ratio.
	 */
	double
	get_llr_threshold() {
		return llr_threshold;
	}


	/// Handle a calculated Z for one strand of an occurrence of a W-mer
    GCC_DIAG_OFF(GCC_UNINIT_WARNING);
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
		if( using_sparse_Z ) {
			if( RevComp ) {
				Z_sparse[ global_pos ].ref().second = Z_strand;
			} else {
				Z_sparse[ global_pos ].ref().first = Z_strand;
			}
		} else {
			if( RevComp ) {
				Z_dense[ global_pos ].second = Z_strand;
			} else {
				Z_dense[ global_pos ].first = Z_strand;
			}
		}

		// check if Z is large enough such that it might have to be normalised later.
		if( Z_strand > large_Z_threshold ) {
			++num_Z_large;
			large_Z_indices.push_back( global_pos );
		}

		// update log likelihood
		LL +=
			safe_x_log_y(
				Z_strand,
				std::log( 1 - g_n ) + log_lambda + log_p_X_n_theta_BS
			)
			+ safe_x_log_y(
				1. - Z_strand,
				std::log( ( 1 - g_n ) * ( 1 - lambda ) + g_n ) + log_p_background
			);
#ifdef MS_DO_DEBUG_CHECKS
		const bool is_finite = ::boost::math::isfinite( LL );
#endif //MS_DO_DEBUG_CHECKS
		MS_CHECK_IS_FINITE( LL, "LL" );
	}
    GCC_DIAG_ON(GCC_UNINIT_WARNING);


	/// Handle a calculated Z for both strands of an occurrence of a W-mer
	inline
	void
	handle_W_mer_all_occurrences(
		top_down_it it,
		size_t num_occurrences,
		double Z_positive_total,
		double Z_negative_total
	) {
		expected_sites += Z_positive_total + Z_negative_total;
		num_Z_non_zero += num_occurrences;

		// update the model
		const prefix_t pre = prefix( representative( it ), base_t::W );
		if( Z_positive_total ) { // first update the positive strand
			base_t::bs_model_t::update_new_model(
				new_bs_model,
				pre,
				Z_positive_total
			);
		}
		if( Z_negative_total ) { // then update the negative strand
			base_t::bs_model_t::update_new_model(
				new_bs_model,
				make_rev_comp_range( pre ),
				Z_negative_total
			);
		}
	}


	/// Check our Z pair has reasonable values.
	static
	inline
	bool
	check_Z( Z_pair_t z ) {
		return
			( ! z.first || is_between_0_and_1( *z.first ) )
			&& ( ! z.second || is_between_0_and_1( *z.second ) )
			;
	}


	/** Adjust a Z. */
	template< typename WMer >
	inline
	double
	adjust_Z( WMer w_mer, double unadjusted_z, double scale ) {
		const double adjusted_z = scale * unadjusted_z;
		const double z_diff = adjusted_z - unadjusted_z;

		// adjust our BS model
		base_t::bs_model_t::update_new_model(
			new_bs_model,
			w_mer,
			z_diff
		);

		// adjust our expected sites
		expected_sites += z_diff;

		// \todo could adjust log likelihood here but would need to recalculate log_p_binding and log_p_background
		//LL -= safe_x_log_y( z_diff, log_lambda_ratio + log_p_binding - log_p_background );

		return adjusted_z;
	}


	/** Normalise a window at global position, j. */
    GCC_DIAG_OFF(GCC_UNINIT_WARNING);
	inline
	bool
	normalise_window_dense( Z_dense_size_t j, Z_dense_size_t max_o ) {
		double local_z = 0.;
		double max_z = 0.; /* find largest z_ij */
		Z_dense_size_t max_p; // and its position
		bool rev_comp; // whether on rev comp or not
		const Z_dense_size_t last_p = std::min( j + base_t::W, max_o );

		for( Z_dense_size_t p = j; p < last_p; p++ ) { /* position */
			const Z_pair_t & Zp = Z_dense[ p ];
			check_Z( Zp );
			MYRRH_ASSERT( check_Z( Zp ) );
			if( Zp.first ) {
				const double z = *Zp.first;
				MYRRH_ASSERT( is_between_0_and_1( z ) );
				local_z += z; /* compute local motif z sum */
				if( z && z > max_z ) { // check if this is the largest z
					max_z = z; /* largest z in window */
					max_p = p; /* position with largest z */
					rev_comp = false;
				}
			}
			if( Zp.second ) {
				const double z = *Zp.second;
				MYRRH_ASSERT( is_between_0_and_1( z ) );
				local_z += z; /* compute local motif z sum */
				if( z && z > max_z ) { // check if this is the largest z
					max_z = z; /* largest z in window */
					max_p = p; /* position with largest z */
					rev_comp = true;
				}
			}
		}
		MYRRH_ASSERT( max_z > 0. ); // make sure we found at least one entry in the window.
		MYRRH_ASSERT( local_z > 0. ); // make sure we found at least one entry in the window.

		/* normalize if necessary; leave largest z in window unchanged */
		if( local_z > 1. ) { /* normalize */
			double scale = (1. - max_z) / (local_z - max_z);
			for( Z_dense_size_t p = j; p < last_p; p++ ) { /* position */
				Z_pair_t & Zp = Z_dense[ p ];
				MYRRH_ASSERT( check_Z( Zp ) );
				if( Zp.first && ( p != max_p || rev_comp ) ) { // don't scale max Z in window
					// reset Z to adjusted value
					Zp.first = adjust_Z( base_t::_data.get_W_mer( base_t::W, p ), *Zp.first, scale );
				}
				if( Zp.second && ( p != max_p || ! rev_comp ) ) { // don't scale max Z in window
					// reset Z to adjusted value
					Zp.second = adjust_Z( make_rev_comp_range( base_t::_data.get_W_mer( base_t::W, p ) ), *Zp.second, scale );
				}
				MYRRH_ASSERT( check_Z( Zp ) );
			} /* position */
			return true;
		} /* normalize */ else {
			return false;
		}
	}


	/** Normalise a window at global position, j. */
	inline
	bool
	normalise_window_sparse( Z_sparse_size_t j, Z_sparse_it j_it, Z_sparse_size_t max_o ) {
		double local_z = 0.;
		double max_z = 0.; /* find largest z_ij */
		Z_sparse_size_t max_p; // and its position
		bool rev_comp; // whether on rev comp or not
		const Z_sparse_size_t last_p = std::min( j + base_t::W, max_o );

		// work out what the local Z is so that we know
		// if the window needs normalising.
		Z_sparse_it end = Z_sparse.end();
		for( Z_sparse_it p_it = j_it; end != p_it && p_it.index() < last_p; ++p_it ) { /* position */
			MYRRH_ASSERT( check_Z( *p_it ) );
			if( (*p_it).first ) {
				const double z = *((*p_it).first);
				MYRRH_ASSERT( is_between_0_and_1( z ) );
				local_z += z; /* compute local motif z sum */
				if( z > max_z ) { // check if this is the largest z
					max_z = z; /* largest z in window */
					max_p = p_it.index(); /* position with largest z */
					rev_comp = false;
				}
			}
			if( (*p_it).second ) {
				const double z = *((*p_it).second);
				MYRRH_ASSERT( is_between_0_and_1( z ) );
				local_z += z; /* compute local motif z sum */
				if( z > max_z ) { // check if this is the largest z
					max_z = z; /* largest z in window */
					max_p = p_it.index(); /* position with largest z */
					rev_comp = true;
				}
			}
		}
		MYRRH_ASSERT( max_z > 0. ); // make sure we found at least one entry in the window.
		MYRRH_ASSERT( local_z > 0. ); // make sure we found at least one entry in the window.

		/* normalize if necessary; leave largest z in window unchanged */
		if( local_z > 1. ) { /* normalize */
			double scale = (1. - max_z) / (local_z - max_z);
			MYRRH_ASSERT( scale < 1. );
			for( Z_sparse_it p_it = j_it; end != p_it && p_it.index() < last_p; ++p_it ) { /* position */
				MYRRH_ASSERT( check_Z( *p_it ) );
				if( (*p_it).first && ( p_it.index() != max_p || rev_comp ) ) {
					// reset Z to adjusted value
				    // get non-const reference
				    Z_sparse_vec::value_type & v = const_cast< Z_sparse_vec::value_type & >( *p_it );
					*(v.first) = adjust_Z( base_t::_data.get_W_mer( base_t::W, p_it.index() ), *((*p_it).first), scale );
				}
				if( (*p_it).second && ( p_it.index() != max_p || ! rev_comp ) ) {
					// reset Z to adjusted value
                    // get non-const reference
                    Z_sparse_vec::value_type & v = const_cast< Z_sparse_vec::value_type & >( *p_it );
					*(v.second) = adjust_Z( make_rev_comp_range( base_t::_data.get_W_mer( base_t::W, p_it.index() ) ), *((*p_it).second), scale );
				}
				MYRRH_ASSERT( check_Z( *p_it ) );
			} /* position */
			return true; // we did normalise
		} else {
			return false; // we did not normalise
		} /* normalize */
	}
    GCC_DIAG_ON(GCC_UNINIT_WARNING);



	/**
	 * Normalise the Z in the window of width W starting at j such that the sum is no more than 1.
	 *
	 * \return next_j
	 *
	 * See the function smooth in the MEME source file tcm.c.
	 */
	Z_sparse_size_t
	normalise_sparse_window_starting_at( Z_sparse_size_t j, Z_sparse_it j_it, Z_sparse_size_t max_o ) {

		// Rewind iterator to the first point at j or after.
		MYRRH_ASSERT( j <= j_it.index() );
        Z_sparse_it begin = Z_sparse.begin();
		while( true ) {
			if( begin == j_it ) { // are we at the beginning?
				break;
			}
			--j_it;
			if( j > j_it.index() ) { // have we gone too far?
				++j_it;
				MYRRH_ASSERT( j_it.index() >= j );
				break;
			}
		}
		MYRRH_ASSERT( j_it.index() >= j );

		// normalise window
		if( normalise_window_sparse( j, j_it, max_o ) ) {
			++num_Z_normalised;
		}

		// return next_j
		return j + base_t::W;
	}



	/**
	 * Normalise the Z in the window of width W starting at j such that the sum is no more than 1.
	 *
	 * \return next_j
	 *
	 * See the function smooth in the MEME source file tcm.c.
	 */
	Z_dense_size_t
	normalise_dense_window_starting_at( Z_dense_size_t j, Z_dense_size_t max_o ) {

		// normalise window
		if( normalise_window_dense( j, max_o ) ) {
			++num_Z_normalised;
		}

		// return next_j
		return j + base_t::W;
	}



	/**
	 * Normalise the Z such that no consecutive W of them sum to more than 1.
	 *
	 * See the function smooth in the MEME source file tcm.c. This is in comments below.
	 */
	void
	normalise_Z() {
		const size_t max_o = Z_sparse.size();  // largest possible offset

		/* normalize adjacent windows of length w, then shift and repeat */
		for( size_t ioff = 0; ioff < std::min( base_t::W, max_o ); ++ioff ) {

			size_t next_j = ioff;

			// Iterate over those Z that are large
			std::sort( large_Z_indices.begin(), large_Z_indices.end() );
			for( index_vec::iterator large_it = large_Z_indices.begin(); large_Z_indices.end() != large_it; ++large_it ) {
				const size_t index = *large_it;

				//if we have already examined a window that covers this index then ignore it
				if( next_j > index ) {
					continue;
				}

				if( using_sparse_Z ) {
					const Z_sparse_it i = Z_sparse.find( index );
					MYRRH_ASSERT( i != i().end () );
					MYRRH_ASSERT( i.index() == index ); // must be in our sparse vector.

					// work out where the beginning of the the window that covers the index should be.
					const Z_sparse_size_t j = index - ((Z_sparse_diff_t(index) - Z_sparse_diff_t(ioff)) % base_t::W);
					MYRRH_ASSERT( next_j <= j );

					next_j = normalise_sparse_window_starting_at( j, i, max_o );
				} else {
					const Z_dense_size_t j = index - ((Z_dense_diff_t(index) - Z_dense_diff_t(ioff)) % base_t::W);
					next_j = normalise_dense_window_starting_at( j, max_o );
				}
			} /* adjacent windows */
		} /* window start */
	}

//		/**
//		 * Normalise the Z such that no consecutive W of them sum to more than 1.
//		 *
//		 * See the function smooth in the MEME source file tcm.c.
//		 */
//		void
//		normalise_Z() {
//			const Z_sparse_size_t max_o = Z.size();  // largest possible offset
//
//			/* normalize adjacent windows of length w, then shift and repeat */
//			for( Z_sparse_size_t ioff = 0; ioff < std::min( base_t::W, max_o ); ioff += 2 ) {
//				for( Z_sparse_size_t j = ioff; j < max_o; j += base_t::W ) { /* adjacent windows */
//					normalise_window( j, max_o );
//				} /* adjacent windows */
//			} /* window start */
//		}
};







} // namespace steme

#endif /* STEME_JR_13AUG2011_EM_H_ */
