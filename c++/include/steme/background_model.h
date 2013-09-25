/** Copyright John Reid 2011, 2012
 *
 * \file
 * \brief Defines the background model type for STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_BACKGROUND_MODEL_H_
#define STEME_JR_13AUG2011_BACKGROUND_MODEL_H_

#include <steme/seqan_types.h>
#include <steme/markov.h>

#include <boost/function_output_iterator.hpp>


namespace steme {


/**
 * Holds likelihoods for every sequence in the text.
 */
typedef std::vector< double >              likelihoods_vec_t;        ///< Holds likelihoods for every character in a sequence.
typedef std::vector< likelihoods_vec_t >   likelihoods_vec_vec_t;    ///< Holds likelihoods for every character in several sequences.
typedef boost::shared_ptr< likelihoods_vec_vec_t > likelihoods_vec_vec_ptr; ///< Pointer to likelihoods.


/// Log likelihood of a W-mer
inline
double
W_mer_log_likelihood( const likelihoods_vec_t & log_likelihoods, size_t offset, size_t W ) {
    if( log_likelihoods.size() <= offset + W - 1 ) {
        throw std::invalid_argument(
            MS_MAKE_STRING(
                "Offset out of range: "
                << log_likelihoods.size() << " <= "
                << offset << " + " << W << " - 1"
            )
        );
    }
    return log_likelihoods[ offset + W - 1 ] - ( offset ? log_likelihoods[ offset - 1 ] : 0. );
}






/**
 * Stores minimum likelihoods for each node in a suffix tree up to depth W.
 */
struct wmer_likelihoods
: boost::noncopyable // don't copy large vector around
{
    typedef boost::shared_ptr< wmer_likelihoods >  ptr;  ///< Pointer.
    typedef seqan::String< double >    likelihoods_vec; ///< Type to hold the likelihoods and their bounds.

    likelihoods_vec            min_log;        ///< Minimum log likelihood of all W-mers under this node.
    double                    LL;             ///< Part of the expected log likelihood under assumption all Z[n] are 0. Used to speed up calculation of LL.
    size_t                     W;              ///< Width of the W-mers we are evaluating.

    /// Constructor.
    inline
    wmer_likelihoods( size_t W )
    : LL( 0. )
    , W( W )
    {
    }

    /// Return the minimum log likelihood of any W-mer that starts with the prefix represented by this iterator.
    template< typename TopDownIt >
    double
    get_min_log_likelihood( TopDownIt it ) const {
        return seqan::getProperty( min_log, value( it ) );
    }
};



/**
 * A background model for a particular motif width, W. This model stores the log likelihoods for each W-mer so
 * is not suitable when there are many W-mers, e.g. when the indexed text is a genome. This class also optionally
 * stores the log likelihood for each occurrence of a W-mer.
 */
struct bg_model
{
    typedef boost::shared_ptr< bg_model >  ptr;  ///< Pointer.

    /// An empty partial evaluation
    struct partial_evaluation_t { };

    wmer_likelihoods               wmer_LLs;       ///< Log likelihoods of the W-mers.
    const likelihoods_vec_vec_t * base_LLs;       ///< Log likelihoods per-base.
    zero_order_frequencies         freqs;          ///< The 0-order frequencies.

    bg_model( size_t W, const zero_order_frequencies & freqs, const likelihoods_vec_vec_t * base_LLs = 0 )
    : wmer_LLs( W )
    , base_LLs( base_LLs )
    , freqs( freqs )
    {
    }


    /// Do we have base-pair resolution log-likelihoods?
    inline
    bool
    have_base_pair_LLs() const {
        return base_LLs;
    }


    /// Return a null partial evaluation.
    static
    partial_evaluation_t
    empty_partial_evaluation() {
        return partial_evaluation_t();
    }


    /// Do nothing as we have no partial evaluation.
    template< typename Prefix >
    partial_evaluation_t
    update_partial_eval( Prefix, partial_evaluation_t, size_t ) {
        return partial_evaluation_t();
    }


    /// Return the minimum log likelihood of any W-mer that starts with the prefix represented by this iterator.
    template< typename TopDownIt >
    double
    get_min_log_likelihood( TopDownIt it, partial_evaluation_t, size_t, size_t ) const {
        return wmer_LLs.get_min_log_likelihood( it );
    }


    /// Return the log likelihood of the W-mer.
    inline
    double
    wmer_log_likelihood( size_t seq, size_t offset ) const {
        if( ! have_base_pair_LLs() ) {
            throw std::logic_error( "Do not have base-pair resolution log likelihoods.");
        }
        if( base_LLs->size() <= seq ) {
            throw std::invalid_argument( "Sequence out of range. " );
        }
        return W_mer_log_likelihood( ( *base_LLs )[ seq ], offset, wmer_LLs.W );
    }


    /// Return the log likelihood of the W-mer if we have it otherwise return the log likelihood of the W-mers indexed by the iterator.
    template< typename TopDownIt >
    double
    wmer_log_likelihood( TopDownIt it, size_t seq, size_t offset ) const {
        // if we have a per-base likelihood use that
        if( have_base_pair_LLs() ) {
            return wmer_log_likelihood( seq, offset );
        } else {
            // otherwise return the smallest likelihood for that W-mer across entire data
            return get_min_log_likelihood( it );
        }
    }
};



/**
 * Calculates likelihoods of W-mers in the data given a model. Uses CRTP C++ idiom.
 * The Derived class is the model.
 */
template< typename Derived >
struct likelihood_calculator
{
    wmer_likelihoods & likelihoods;

    likelihood_calculator( wmer_likelihoods & likelihoods )
    : likelihoods( likelihoods )
    {
    }


	/** Does the calculation. */
	template< typename Data >
	void
	calculate( Data & _data ) {
		using namespace seqan;
		resizeVertexMap( _data.index, likelihoods.min_log );
		visit_node( _data, static_cast< Derived * >( this )->index_root( _data ), 0 );
	}



    /** The iterator used for the descent. Can be overridden in derived class. */
    template< typename Data >
    typename Data::top_down_it
    index_root( Data & _data ) {
        return _data.index_root();
    }


	/**
	 * Calculates likelihood of W-mers under the model.
	 *
	 * \return Minimum likelihood bound for nodes below.
	 */
	template<
		typename Data,
		typename It
	>
	double
	visit_node( Data & _data, It it, size_t last_length ) {
		using namespace seqan;

		MS_DEBUG_STRING_VAR( w_mer, representative( it ) );

		// index into our bounds.
		typename VertexDescriptor<
			typename Container< It >::Type
		>::Type v = value( it );

		double min_ll = 0.; // Return a high bound on the log likelihood if we don't find a lower one.

		// how far down the tree are we?
		const size_t w = std::min( repLength( it ), likelihoods.W );

		// only deal with nodes that do not have 'N's
		if( ! _data.contains_n( likelihoods.W, it ) ) {

	        static_cast< Derived * >( this )->entered_node( it, last_length );

            // did we get to the bottom?
            if( likelihoods.W == w ) {

                // yes, we are at bottom
                // so calculate the log likelihood
                double W_mer_min_ll, total_ll;
                boost::tie( W_mer_min_ll, total_ll ) = static_cast< Derived * >( this )->log_likelihoods_for_w_mer( it );

                // update our calculation of part of the LL
                likelihoods.LL += total_ll;
                min_ll = std::min( min_ll, W_mer_min_ll );

            } else {

                // continue : so go down and across
                if( goDown( it ) ) {
                    min_ll = visit_node( _data, it, w );
                    while( goRight( it ) ) {
                        min_ll = std::min( min_ll, visit_node( _data, it, w ) );
                    }
                }
            }

            // store the lowest log likelihood
            MYRRH_ASSERT( min_ll <= 0. ); // perhaps not for empty tree...
            seqan::assignProperty( likelihoods.min_log, v, min_ll );

            static_cast< Derived * >( this )->left_node( it, last_length );
        }

		return min_ll;
	}
};


/**
 * Partition a sequence into segments based on whether a predicate is true or not.
 */
template< typename Iterator, typename Predicate, typename TrueCallback, typename FalseCallback >
void
segmentize_sequence(
	Iterator begin,
	Iterator end,
	Predicate & predicate,
	TrueCallback true_segment_cb,
	FalseCallback false_segment_cb
) {
	Iterator true_segment_start = end;
	Iterator false_segment_start = end;
	// for each character in the sequence.
	Iterator i = begin;
	for( ; end != i; ++i ) {
		// is the predicate true for this character?
		if( ! predicate( i ) ) {
			// yes we have an false value: did we complete a true segment?
			if( end != true_segment_start ) {
				// we have completed a true segment
				MYRRH_ASSERT( end == false_segment_start );
				true_segment_cb( true_segment_start, i );
				true_segment_start = end; // reset beginning of true segment
			}
			// is this the beginning of an true segment?
			if( end == false_segment_start ) {
				false_segment_start = i; // remember start of false segment
			}
		} else {
			// no we have an true value: did we complete a false segment?
			if( end != false_segment_start ) {
				// we have completed a false segment
				MYRRH_ASSERT( end == true_segment_start );
				false_segment_cb( false_segment_start, i );
				false_segment_start = end; // reset beginning of false segment
			}
			// did we start a true segment?
			if( end == true_segment_start ) {
				// yes we did start a true segment
				true_segment_start = i; // remember start of true segment
			}
		}
	}
	// handle last segment
	if( end != true_segment_start ) {
		MYRRH_ASSERT( end == false_segment_start );
		true_segment_cb( true_segment_start, i );
	} else {
		if( begin != end ) { // could be empty input sequence
			MYRRH_ASSERT( end != false_segment_start );
			false_segment_cb( false_segment_start, i );
		}
	}
}

/// Is the value of the iterator an unknown value?
template< typename Iterator >
bool
is_unknown( Iterator i ) {
	return seqan::value( i ) == seqan::unknownValue< typename seqan::Value< Iterator >::Type >();
}

/**
 * Partition a sequence into segments based on whether the bases are known or not.
 */
template< typename Seq, typename KnownCallback, typename UnknownCallback >
void
segmentize_by_known_bases(
	const Seq & seq,
	KnownCallback known_segment_cb,
	UnknownCallback unknown_segment_cb
) {
	typedef typename seqan::Iterator< const Seq >::Type it_t;
	segmentize_sequence(
		seqan::begin( seq ),
		seqan::end( seq ),
		is_unknown< it_t >,
		unknown_segment_cb,
		known_segment_cb
	);
}

/// Calculate the likelihoods of the known segments
template< typename CompleteMarkovModel, typename Iterator >
void
calculate_known_segment_likelihoods(
	const CompleteMarkovModel & mm,
	likelihoods_vec_t & likelihoods,
	Iterator segment_start,
	Iterator segment_end
) {
	mm.evaluate(
		segment_start,
		segment_end,
		boost::make_function_output_iterator( make_push_back_cumulative( likelihoods ) )
	);
}

/// Calculate the likelihoods of the unknown segments.
template< typename Iterator >
void
calculate_unknown_segment_likelihoods(
	likelihoods_vec_t & likelihoods,
	Iterator segment_start,
	Iterator segment_end
) {
	// deal with unknown values
	push_back_cumulative< likelihoods_vec_t > adder( likelihoods );
	for( ; segment_end != segment_start; ++segment_start ) {
		adder( log_quarter );
	}
}

/// Calculate the likelihoods for a sequence
template< typename CompleteMarkovModel, typename Seq  >
void
calculate_likelihoods( const CompleteMarkovModel & mm, const Seq & seq, likelihoods_vec_t & likelihoods ) {
	using boost::lambda::bind;
	typedef typename seqan::Iterator< const Seq >::Type it_t;
	likelihoods.clear();
	segmentize_by_known_bases(
		seq,
		boost::lambda::bind(
			&calculate_known_segment_likelihoods< CompleteMarkovModel, it_t >,
			boost::ref( mm ),
			boost::ref( likelihoods ),
			boost::lambda::_1,
			boost::lambda::_2
		),
		boost::lambda::bind(
			&calculate_unknown_segment_likelihoods< it_t >,
			boost::ref( likelihoods ),
			boost::lambda::_1,
			boost::lambda::_2
		)
	);
	// make sure we calculated the right number of likelihoods
	MYRRH_ASSERT( likelihoods.size() == seqan::length( seq ) );
}

/// Calculate the likelihoods from the sequences
template< typename CompleteMarkovModel, typename Text >
void
calculate_likelihoods( const CompleteMarkovModel & mm, const Text & text, likelihoods_vec_vec_t & likelihoods ) {
	likelihoods.clear();
	for( size_t i = 0; seqan::length( text ) != i; ++i ) {
		likelihoods.push_back( likelihoods_vec_t() );
		calculate_likelihoods( mm, seqan::getValue( text, i ), likelihoods[i] );
	}
}

/// Adds the known segment counts to the Markov model.
template< typename CompleteMarkovModel >
struct add_known_segment_counts {
	CompleteMarkovModel & mm;
	add_known_segment_counts( CompleteMarkovModel & mm ) : mm( mm ) { } ///< Constructor.

	/// Add the counts
	template< typename Iterator >
	void operator()( Iterator begin, Iterator end ) {
		add_counts_to_complete( mm, begin, end );
	}
};

/// Add the counts from the unknown segments
struct add_unknown_segment_counts {
	/// Add the counts
	template< typename Iterator >
	void operator()( Iterator begin, Iterator end ) {
	}
};

template< typename CompleteMarkovModel, typename Text >
void
add_counts_from_seqs_to_complete( CompleteMarkovModel & model, const Text & text ) {
	// add the counts from the sequences
	for( size_t i = 0; seqan::length( text ) != i; ++i ) {
		segmentize_by_known_bases(
			seqan::getValue( text, i ),
			add_known_segment_counts< CompleteMarkovModel >( model ),
			add_unknown_segment_counts()
		);
	}
}

/// Build a Markov model of the sequences and return the 0-order frequencies
template< typename CompleteMarkovModel, typename Text >
zero_order_frequencies
build_model_from_text(
	CompleteMarkovModel & model,
	const Text & text,
	double pseudo_count
) {
	add_counts_from_seqs_to_complete( model, text );

	zero_order_frequencies freqs(
		model.template get_lower_order_model< 0 >().mm.storage.begin(),
		model.template get_lower_order_model< 0 >().mm.storage.end()
	);

	// add the pseudo-counts, normalise and convert to the correct scale.
	add_pseudo_counts( model, pseudo_count );
	normalise_counts( model );
	convert_to_scale( model );

	return freqs;
}



namespace detail {

/// Add the counts from the index to the markov model
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct add_counts_from_index_descender
: descend_to_depth< add_counts_from_index_descender< Order, AlphabetSize, Float, Scale > >
{
	typedef markov_model< Order, AlphabetSize, Float, Scale > mm_t; ///< The Markov model type.

	mm_t & mm; ///< The Markov model.

	/// Constructor
	add_counts_from_index_descender( mm_t & mm )
	: descend_to_depth< add_counts_from_index_descender >( Order + 1 )
	, mm( mm )
	{ }

	/// add the counts for the W-mers at this iterator
	template< typename TopDownIt >
	void
	handle_node( TopDownIt it ) {

		using namespace seqan;

		const int first_unknown = find_first_unknown( prefix( representative( it ), Order + 1 ) );

		// if no unknown values in prefix then
		if( -1 == first_unknown ) {
			detail::access_multi_array_element( boost::type< Float & >(), mm.storage, begin( representative( it ) ) )
				+= countOccurrences( it )
				;
		}
	}
};

/// Add counts to the Markov model from the index.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale, typename Index >
void
add_counts_from_index( markov_model< Order, AlphabetSize, Float, Scale > & mm, Index & index ) {
	add_counts_from_index_descender< Order, AlphabetSize, Float, Scale >( mm ).descend_index( index );
}

/// Add counts from an index to a complete Markov model.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct
add_counts_from_index_to_complete {
	template< typename Index >
	void
	operator()( complete_markov_model< Order, AlphabetSize, Float, Scale > & model, Index & index ) const {
		// add counts for this order
		add_counts_from_index( model.mm, index );
		// add counts for lower orders
		add_counts_from_index_to_complete< Order-1, AlphabetSize, Float, Scale >()( model.lower_orders, index );
	}
};

/// Specialisation for 0-order
template< size_t AlphabetSize, typename Float, typename Scale >
struct
add_counts_from_index_to_complete< 0, AlphabetSize, Float, Scale > {
	template< typename Index >
	void
	operator()( complete_markov_model< 0, AlphabetSize, Float, Scale > & model, Index & index ) const {
		// add counts for this order
		add_counts_from_index( model.mm, index );
	}
};


} //namespace detail


/// Add counts from sequence to complete model.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale, typename Index >
void
add_counts_from_index_to_complete( complete_markov_model< Order, AlphabetSize, Float, Scale > & model, Index & index ) {
	detail::add_counts_from_index_to_complete< Order, AlphabetSize, Float, Scale >()( model, index );
}




/// Build a Markov model of the sequences in the index and return the 0-order frequencies
template<
	typename CompleteMarkovModel,
	typename Index
>
zero_order_frequencies
build_model_from_index(
	CompleteMarkovModel & model,
	Index & index,
	double pseudo_count
) {
	// add the counts from the index
	add_counts_from_index_to_complete( model, index );

	// get the zero-order frequencies
	zero_order_frequencies freqs(
		model.template get_lower_order_model< 0 >().mm.storage.begin(),
		model.template get_lower_order_model< 0 >().mm.storage.end()
	);

	// add the pseudo-counts, normalise and convert to the correct scale.
	add_pseudo_counts( model, pseudo_count );
	normalise_counts( model );
	convert_to_scale( model );

	return freqs;
}



/**
 * Markov 0-order background model.
 */
struct markov_0_order_likelihoods_calculator
: likelihood_calculator< markov_0_order_likelihoods_calculator >
{
    typedef likelihood_calculator< markov_0_order_likelihoods_calculator > base_t;
    typedef boost::shared_ptr< markov_0_order_likelihoods_calculator >  ptr;                ///< A pointer to a Markov0BackgroundModel.

    zero_order_frequencies         freqs;          ///< The 0-order frequencies.

    /** Constructor. */
    template< typename Data >
    markov_0_order_likelihoods_calculator(
        Data & _data,
        const zero_order_frequencies & freqs,
        wmer_likelihoods & likelihoods
    )
    : base_t( likelihoods )
    , freqs( freqs )
    {
        base_t::calculate( _data );
    }


    /**
     * Called by base class
     */
    template< typename It >
    void
    entered_node( It it, size_t last_length ) {
    }

    /**
     * Called by base class
     */
    template< typename It >
    void
    left_node( It it, size_t last_length ) {
    }


    /** Return the minimum and total log likelihoods of the occurrences of this prefix. */
    template< typename It >
    boost::tuple< double, double >
    log_likelihoods_for_w_mer( It i ) {
        using namespace seqan;

        typedef typename Container< It >::Type index_t;
        typedef typename Infix< typename Fibre< index_t, EsaText >::Type const >::Type representative_t;
        typedef typename Prefix< representative_t >::Type prefix_t;
        typedef typename Iterator< prefix_t >::Type prefix_it;
        typedef typename Value< prefix_it >::Type base_t;

        prefix_t p = prefix( representative( i ), likelihoods.W );
        MS_DEBUG_STRING_VAR( w_mer, p );
        double ll = 0.;
        for( prefix_it it = begin( p ); it != end( p ); ++it ) {
            base_t base = value( it );
            if( base == unknownValue< base_t >() ) {
                ll += log_quarter; // treat unknowns as log 1/4
            } else {
                ll += freqs.dist_logs[ ordValue( base ) ];
            }
        }

        return boost::make_tuple( ll, countOccurrences( i ) );
    }
};





/**
 * A background model created from likelihoods for each base in the sequences.
 */
struct base_to_wmer_likelihoods_calculator
: likelihood_calculator< base_to_wmer_likelihoods_calculator >
{
    typedef likelihood_calculator< base_to_wmer_likelihoods_calculator > base_t;
	typedef boost::shared_ptr< base_to_wmer_likelihoods_calculator >  ptr;  ///< A pointer to a BackgroundModelFromLikelihoods.

	const likelihoods_vec_vec_t &        base_likelihoods;     ///< The log likelihoods for every base in the sequences.

	/** Constructor. */
	template< typename Data >
	base_to_wmer_likelihoods_calculator(
		Data & _data,
		const likelihoods_vec_vec_t & base_likelihoods,
        wmer_likelihoods & likelihoods
	)
	: base_t( likelihoods )
	, base_likelihoods( base_likelihoods )
	{
	    base_t::calculate( _data );
	}

    /**
     * Called by base class
     */
    template< typename It >
    void
    entered_node( It it, size_t last_length ) {
    }

    /**
     * Called by base class
     */
    template< typename It >
    void
    left_node( It it, size_t last_length ) {
    }


	/** Likelihoods for W-mer. */
	inline
	double
	w_mer_log_likelihood( size_t seq, size_t offset ) {
	    if( base_likelihoods.size() <= seq ) {
	        throw std::invalid_argument( "Sequence out of range. " );
	    }
	    return W_mer_log_likelihood( base_likelihoods[ seq ], offset, likelihoods.W );
	}


	/** Return the minimum and total log likelihoods of the occurrences of this prefix. */
	template< typename It >
	boost::tuple< double, double >
	log_likelihoods_for_w_mer( It it ) {
		using namespace seqan;

		typedef typename Container< It >::Type index_t;
		typedef typename Infix< typename Fibre< index_t, EsaSA >::Type const >::Type occurrences_t;
		typedef typename Value< occurrences_t >::Type occurrence_t;
		occurrences_t occs = seqan::getOccurrences( it );
		double min_ll = 0.;
		double total_ll = 0.;
		for( unsigned i = 0; seqan::length( occs ) != i; ++i ) {
			const occurrence_t & occ = occs[ i ];
			const size_t seq = getSeqNo( occ );
			const size_t off = getSeqOffset( occ );
			const double likelihood = w_mer_log_likelihood( seq, off );
			total_ll += likelihood;
			min_ll = std::min( likelihood, min_ll );
			MYRRH_ASSERT( min_ll < 0. );
		}
		return boost::make_tuple( min_ll, total_ll );
	}
};




/**
 * Calculates the likelihoods of the W-mers using a Markov model.
 */
template< typename MarkovModel >
struct markov_wmer_likelihoods_calculator
: likelihood_calculator< markov_wmer_likelihoods_calculator< MarkovModel > >
{
    typedef likelihood_calculator< markov_wmer_likelihoods_calculator< MarkovModel > > base_t;
    typedef MarkovModel markov_model_t;

    const markov_model_t &    model;
    std::vector< double >     current_evals;

    template< typename Data >
    markov_wmer_likelihoods_calculator(
        Data & _data,
        const MarkovModel & model,
        wmer_likelihoods & likelihoods
    )
    : base_t( likelihoods )
    , model( model )
    {
        base_t::calculate( _data );
    }


    /**
     * Called by base class
     */
    template< typename It >
    void
    entered_node( It it, size_t last_length ) {
        using namespace seqan;
        MS_DEBUG_STRING_VAR( w_mer, representative( it ) );

        // work out where we are
        const size_t length = std::min( base_t::likelihoods.W, repLength( it ) );

        // resize evaluations
        MYRRH_ASSERT( current_evals.size() >= last_length );
        current_evals.resize( last_length );

        // evaluate part of the representative
        typedef typename Container< It >::Type IndexT;
        typedef typename Infix< typename Fibre< IndexT, EsaText >::Type const>::Type ReprT;
        typedef typename boost::range_iterator< ReprT >::type range_it;
        range_it begin = boost::begin( representative( it ) );
        model.evaluate_in_context(
            begin,
            begin + last_length,
            begin + length,
            boost::make_function_output_iterator( make_push_back_cumulative( current_evals ) )
        );
    }

    /**
     * Called by base class
     */
    template< typename It >
    void
    left_node( It it, size_t last_length ) {
    }

    /**
     * Return the minimum and total log likelihoods of the occurrences of this prefix.
     */
    template< typename It >
    boost::tuple< double, double >
    log_likelihoods_for_w_mer( It it ) {
        using namespace seqan;
        MS_DEBUG_STRING_VAR( w_mer, representative( it ) );
        const double it_LL = current_evals.back();
        return boost::make_tuple( it_LL, countOccurrences( it ) * it_LL );
    }
};



/// Meta-function to choose type.
template< typename Spec = default_spec >
struct background_model_meta {
    typedef bg_model type; ///< Background model type.
};




} // namespace steme

#endif // STEME_JR_13AUG2011_BACKGROUND_MODEL_H_
