/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the binding model type for STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_BINDING_MODEL_H_
#define STEME_JR_13AUG2011_BINDING_MODEL_H_

#include <steme/defs.h>
extern "C" {
#include <steme/meme.h>
}

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/shared_ptr.hpp>




namespace steme {

/** Whether MEME LLR p-value tables have been initialised. */
extern boost::optional< boost::tuple< int, int > > MEME_llr_pv_tables_initialised;




/** A storage for (log) probabilities or weights in a PSSM. */
struct PssmStorage
{
	typedef boost::shared_ptr< PssmStorage >            ptr;                 ///< Pointer to this type.
	typedef boost::multi_array< double, 2 >             matrix;              ///< Our matrix type.
	typedef matrix::size_type                           size_type;           ///< Size type.

	matrix                m;                ///< Storage for probabilities, weights or log probabilities.

	/** Constructor. */
	PssmStorage( size_t W, size_t alphabet_size = 4 )
	: m( boost::extents[W][alphabet_size] )
	{
	}

	/** Copy constructor. */
	PssmStorage( const PssmStorage & other )
	: m( other.m )
	{
	}

	/** Width of PSSM. */
	size_type
	W() const {
		return boost::size( m );
	}

	/** Size of the alphabet. */
	size_type
	alphabet_size() const {
		return m.shape()[1];
	}

	PssmStorage
	rev_comp() const {
		PssmStorage result( W(), alphabet_size() );
		for( size_type w = 0; W() != w; ++w ) {
			for( size_type b = 0; alphabet_size() != b; ++b ) {
				result.m[W()-w-1][alphabet_size()-b-1] = m[w][b];
			}
		}
		return result;
	}
};




/** Apply a functor to the elements of a multi_array. */
template< typename MultiArray, typename Fn >
void
apply_functor( MultiArray & m, Fn f )
{
	BOOST_FOREACH(
		typename MultiArray::element & p,
		boost::make_iterator_range( m.data(), m.data() + m.num_elements() )
	) {
		p = f( p );
	}
}



/**
 * Normalise a Matrix.
 *
 * \return How many counts in its last column.
 */
template< typename Matrix >
double
normalise( Matrix & m )
{
	double total = 0.;
	BOOST_FOREACH( typename Matrix::reference col, m ) {
		MYRRH_ASSERT( 4 == boost::size( col ) );
		total = std::accumulate( col.begin(), col.end(), 0.0 );
		if( total <= 0. ) {
			throw std::logic_error( "No counts to normalise." );
		}
		BOOST_FOREACH( double & p, col ) {
			p /= total;
		}
	}
	return total;
}



/**
 * Evaluate the PSSM on the (partial) sequence (or its reverse complement).
 *
 * @arg pssm_log_probs: The log probabilities of the bases under the PSSM.
 * @arg start_w: Start at this index in the PSSM. This is used when updating old partial evaluations.
 * @arg threshold: If the evaluation goes below this threshold, then stop evaluating.
 *
 * @todo Handle 'N's
 */
template< bool RevComp, typename Seq >
boost::optional< double >
evaluate_pssm(
	Seq & seq,
	const PssmStorage & pssm_log_probs,
	size_t start_w = 0,
	double threshold = -std::numeric_limits< double >::max()
)
{
	MS_DEBUG_STRING_RANGE( dbg_seq, seq );

	typedef PssmStorage::size_type size_type;
	const size_type len = boost::size( seq );
	const size_type W = pssm_log_probs.W();
	const size_type alphabet_size = pssm_log_probs.alphabet_size();
	MYRRH_ASSERT( len <= W );
	MYRRH_ASSERT( start_w <= len );
	double evaluation = 0.0;

	// for each base in the PSSM starting at start_w position
	typename boost::range_iterator< Seq >::type start = boost::begin( seq ) + start_w;
	for( size_type w = start_w; len != w; ++w ) {
		const size_type b = seqan::ordValue( *start );
		MYRRH_ASSERT( 4 >= b ); //check we don't have unknown bases....
		const double log_prob =
			RevComp
				? pssm_log_probs.m[W-w-1][alphabet_size-b-1]
				: pssm_log_probs.m[w][b]
				;
		evaluation += log_prob;
		if( evaluation < threshold ) { // check if we have gone below the threshold
			return boost::optional< double >();
		}
		++start;
	}
	return evaluation;
}



/**
 * I previously had this down as calculating the LLR (as in MEME) of a distribution. However it seems
 * to calculate the information content so I've renamed it.
 *
 * Both the distribution and the background are assumed to be given in the log scale.
 */
template< typename Dist >
double
calculate_column_ic( const Dist & dist_logs, const double * back_logs ) {
	double result = 0.;
	typedef typename boost::range_iterator< const Dist >::type dist_it;
	for( dist_it f_it = boost::begin( dist_logs ); f_it != boost::end( dist_logs ); ++f_it, ++back_logs ) {
		const double log_f = *f_it;
		if( ::boost::math::isfinite( log_f ) ) { //ignore impossible bases
			const double freq = std::exp( log_f );
			const double log_p = *back_logs;
			result += freq * ( log_f - log_p );
		}
	}
	return result;
}



/** A PSSM. */
struct Pssm
{
	typedef boost::shared_ptr< Pssm >        ptr;           ///< A pointer to a PSSM.

	PssmStorage     log_probs;                    ///< The PSSM's log probabilities.
	double          num_samples;                  ///< Number of samples the model is created from.

	/** Constructor. */
	Pssm( PssmStorage log_probs )
		: log_probs( log_probs )
		, num_samples( 0. )
	{
	}

	/** Check the exponentiated log probabilities sum to close to 1. */
	bool
	probabilities_sum_to_1( double tolerance = .01 ) const {
		BOOST_FOREACH( const PssmStorage::matrix::const_reference pssm_col, log_probs.m ) {
			if( ! log_probs_sum_to_1( pssm_col, tolerance ) ) {
				return false;
			}
		}
		return true;
	}


	/**
	 * Calculate the log likelihood ratio of the PSSM.
	 */
	double
	calculate_llr( double * back_dist_logs ) {
		return - num_samples * calculate_information_content( back_dist_logs );
	}


	/**
	 * Calculate the information content of the PSSM.
	 */
	double
	calculate_information_content( double * back_dist_logs ) {
		double ic = 0.;
		BOOST_FOREACH( PssmStorage::matrix::reference pssm_col, log_probs.m ) {
			ic += calculate_column_ic( pssm_col, back_dist_logs );
		}
		return ic;
	}

	/** The reverse complement of the PSSM. */
	Pssm
	rev_comp() const {
		Pssm result( log_probs.rev_comp() );
		result.num_samples = num_samples;
		return result;
	}

};



/**
 * PSSM binding site model. Implements the interface STEME expects.
 */
class PssmBindingSiteModel
{

public:

	typedef boost::shared_ptr< PssmBindingSiteModel >     ptr;                          ///< Pointer to this type.
	typedef std::vector< double >                         vector;                       ///< A vector of doubles.
	typedef boost::shared_ptr< PssmStorage::matrix >      new_model_t;                  ///< Type used for new models.
	typedef double                                        partial_evaluation_t;         ///< The type of a partial evaluation under the model.


public:

	Pssm            pssm;                         ///< The PSSM.
	Pssm            rev_comp_pssm;                ///< A reverse-complemented copy of the PSSM for efficiency.
	double          seed_pseudo_counts;           ///< Pseudo-counts used when seeding the model.
	vector          best_suffixes;                ///< Store probabilities of best suffixes.
	vector          best_prefixes;                ///< Store probabilities of best prefixes.


public:

	/** Constructor. */
	PssmBindingSiteModel( Pssm pssm, double seed_pseudo_counts=.05 )
		: pssm( pssm )
		, rev_comp_pssm( pssm.rev_comp() )
		, seed_pseudo_counts( seed_pseudo_counts )
		, best_suffixes( pssm.log_probs.W() + 1 )
		, best_prefixes( pssm.log_probs.W() + 1 )
	{
		recalculate();
	}

	/** Seed the model with the W-mer. */
	template< typename Prefix >
	void
	seed( Prefix prefix, bool use_pseudo_count = true ) {
		if( seqan::length( prefix ) != pssm.log_probs.W() ) { //check seed is correct length
			throw std::logic_error(
				MYRRH_MAKE_STRING(
					"Seed not correct length. Seed length="<<seqan::length( prefix )<<"; model length="<<pssm.log_probs.W()
				)
			);
		}
		typename boost::range_iterator< Prefix >::type start = boost::begin( prefix );
		for( size_t w = 0; pssm.log_probs.W() != w; ++w, ++start ) {
			std::fill( pssm.log_probs.m[w].begin(), pssm.log_probs.m[w].end(), use_pseudo_count ? seed_pseudo_counts : 0. );
			const size_t b = seqan::ordValue( *start );
			if( 4 < b ) {
				throw std::logic_error( "Strange base in seed" );
			}
			if( 4 != b ) { //ignore 'N's
				pssm.log_probs.m[w][b] += 1.;
			}
		}
		normalise( pssm.log_probs.m );
		double( * log_fn )( double ) = &std::log;
		apply_functor( pssm.log_probs.m, log_fn );
		MYRRH_ASSERT( pssm.probabilities_sum_to_1() );
		recalculate();
	}


	/** Takes a previous partial evaluation and updates it. Thresholded version so that if the evaluation
	 * of the entire W-mer is going to be below the threshold an empty optional is returned. */
	template<
		bool RevComp,
		typename Prefix
	>
	boost::optional< partial_evaluation_t >
	update_partial_eval_thresholded(
		Prefix prefix,
		partial_evaluation_t previous_eval,
		size_t previous_w,
		double threshold
	) const	{
		MS_DEBUG_STRING_RANGE( dbg_seq, prefix );
		const double remainder = best_remainder< RevComp >( seqan::length( prefix ) );
		const double new_threshold = threshold - previous_eval - remainder;
		const boost::optional< double > eval_update =  evaluate_pssm< false >(
			prefix,
			RevComp ? rev_comp_pssm.log_probs : pssm.log_probs,
			previous_w,
			new_threshold
		);
		return eval_update
			? *eval_update + previous_eval
			: boost::optional< partial_evaluation_t >()
			;
	}


	/** Returns a bound on the evaluation. */
	template<
		bool RevComp,
		typename Prefix
	>
	partial_evaluation_t
	bound_evaluation( const Prefix & prefix ) {
		MS_DEBUG_STRING_RANGE( dbg_seq, prefix );
		const boost::optional< double > prefix_eval = evaluate_pssm< false >(
			prefix,
			RevComp ? rev_comp_pssm.log_probs : pssm.log_probs
		);
		MYRRH_ASSERT( prefix_eval ); // we didn't give it a threshold so we must have an evaluation
		return *prefix_eval + best_remainder< RevComp >( seqan::length( prefix ) );
	}


	/** Returns an empty evaluation. */
	static
	inline
	partial_evaluation_t
	empty_partial_evaluation() {
		return 0.;
	}


	/// Takes an evaluation and returns the log likelihood.
	static
	inline
	double
	eval_to_log( partial_evaluation_t eval ) {
		return eval;
	}


	/** Create a new model. */
	new_model_t
	create_new_model( bool use_pseudo_counts = true ) const {
		new_model_t new_model(
			new PssmStorage::matrix( boost::extents[ pssm.log_probs.W() ][ pssm.log_probs.alphabet_size() ] )
		);
		std::fill( new_model->data(), new_model->data() + new_model->num_elements(), use_pseudo_counts ? seed_pseudo_counts : 0. );
		return new_model;
	}


	/** Update the model using the newly learnt model. */
	void
	update_from_new( new_model_t new_model ) {
		pssm.log_probs.m = *new_model;
		pssm.num_samples = normalise( pssm.log_probs.m );
		double( * log_fn )( double ) = &std::log;
		apply_functor( pssm.log_probs.m, log_fn );
		MYRRH_ASSERT( pssm.probabilities_sum_to_1() );
		recalculate();
	}


	/**
	 * Update the new model using the sequence and its evaluation.
	 */
	template< typename Seq >
	static
	void
	update_new_model( new_model_t new_model, const Seq & seq, double weight ) {
		MS_DEBUG_STRING_RANGE( dbg_seq, seq );

		MYRRH_ASSERT( boost::size( seq ) == boost::size( *new_model ) );
		typename boost::range_iterator< const Seq >::type start = boost::begin( seq );
		const size_t len = boost::size( seq );
		for( size_t w = 0; len != w; ++w ) {
			const size_t x = seqan::ordValue( *start );
			(*new_model)[ w ][ x ] += weight;
			++start;
		}
	}


    /** Recalculate best prefixes and suffixes. */
    void
    recalculate() {
        MYRRH_ASSERT( boost::size( best_suffixes ) == size_t( pssm.log_probs.W() + 1 ) );
        best_suffixes[ pssm.log_probs.W() ] = 0.;
        calculate_best_fixes( pssm.log_probs.m.rbegin(), pssm.log_probs.m.rend(), best_suffixes.rbegin() + 1 );

        MYRRH_ASSERT( boost::size( best_prefixes ) == size_t( pssm.log_probs.W() + 1 ) );
        best_prefixes[ pssm.log_probs.W() ] = 0.;
        calculate_best_fixes( pssm.log_probs.m.begin(), pssm.log_probs.m.end(), best_prefixes.rbegin() + 1 );

        rev_comp_pssm = pssm.rev_comp();
    }


protected:
	/** Return the best likelihood for the rest of the model. */
	template< bool RevComp >
	inline
	double
	best_remainder( size_t w ) const {
		return RevComp ? best_prefixes[w] : best_suffixes[w];
	}


	/** Calculate log probabilities of best prefixes or suffixes. */
	template< typename It, typename OutIt >
	static
	void
	calculate_best_fixes( It begin, It end, OutIt out_it ) {
		double best = 0.;
		while( begin != end ) {
			const double maximum = *std::max_element( (*begin).begin(), (*begin).end() );
			best += maximum;
			*out_it = best;
			++out_it;
			++begin;
		}
	}


};


/// Meta-function to choose type.
template< typename Spec = default_spec >
struct binding_model_meta {
	typedef PssmBindingSiteModel type; ///< Binding site model.
};



/// Evaluate a W-mer
template<
    bool RevComp,
    typename BSModel,
    typename Wmer
>
double
w_mer_binding_log_likelihood( const BSModel & bs_model, const Wmer & w_mer )
{
    // update for W-mer
     boost::optional< typename BSModel::partial_evaluation_t > eval = bs_model.template update_partial_eval_thresholded< RevComp >(
        w_mer,
        BSModel::empty_partial_evaluation(),
        0,
        -std::numeric_limits< double >::max()
    );
    BOOST_ASSERT( eval );

    return BSModel::eval_to_log( *eval );
}




}

#endif // STEME_JR_13AUG2011_BINDING_MODEL_H_
