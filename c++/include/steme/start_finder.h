/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the start finder type for STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_START_FINDER_H_
#define STEME_JR_13AUG2011_START_FINDER_H_

#include <steme/find_best_w_mers.h>
#include <steme/find_instances.h>
#include <steme/significance.h>

#include <boost/function.hpp>
#include <boost/range/algorithm/min_element.hpp>

namespace steme {



typedef unsigned long long              seq_ordinal_t;     ///< Holds the result of a sequence ordinal calculation.
typedef std::set< seq_ordinal_t >       seq_ordinal_seq;   ///< A set of sequence ordinals.


/** Checks to see if the ordinal type is large enough to hold ordinals for all the sequences of length w. */
template< size_t alphabet_size, typename ordinal_t >
bool
ordinal_type_large_enough( size_t w ) {
	return w * std::log( alphabet_size ) < std::log( std::numeric_limits< ordinal_t >().max() );
}


/**
 * Calculate unique value for a sequence (shared with reverse complement).
 */
template< size_t alphabet_size, typename SeqRange >
seq_ordinal_t
calculate_seq_ordinal_one_way( const SeqRange & seq ) {
	seq_ordinal_t ordinal = 0;
	seq_ordinal_t power = 1;
	BOOST_FOREACH( typename boost::range_value< SeqRange >::type c, seq ) {
		const seq_ordinal_t ord = seqan::ordValue( c );
		MYRRH_ASSERT( ord < alphabet_size );
		ordinal += ord * power;
		power *= alphabet_size;
	}
	return ordinal;
}


/**
 * Edit distance between 2 strings.
 */
template<
    typename Str1,
    typename Str2
>
unsigned
edit_distance( const Str1 & s1, const Str2 & s2 ) {
    BOOST_ASSERT( boost::size( s1 ) == boost::size( s2 ) );
    unsigned distance = 0;
    for( unsigned i = 0; boost::size( s1 ) != i; ++i ) {
        if( s1[ i ] != s2[ i ] ) {
            ++distance;
        }
    }
    return distance;
}


/**
 * Edit distance between 2 strings including possible shifts.
 */
template<
    typename Str1,
    typename Str2
>
unsigned
shifted_edit_distance( const Str1 & s1, const Str2 & s2 ) {
    BOOST_ASSERT( boost::size( s1 ) == boost::size( s2 ) );
    unsigned distance = edit_distance( s1, s2 );
    for( unsigned i = 1; boost::size( s1 ) > i && distance > i; ++i ) {
        distance = std::min(
            i + edit_distance(
                boost::make_iterator_range( s1.begin() + i, s1.end() ),
                boost::make_iterator_range( s2.begin(), s2.end() - i )
            ),
            distance
        );
        distance = std::min(
            i + edit_distance(
                boost::make_iterator_range( s1.begin(), s1.end() - i ),
                boost::make_iterator_range( s2.begin() + i, s2.end() )
            ),
            distance
        );
    }
    return distance;
}




/**
 * Possibly reverse-complemented edit distance between 2 strings including possible shifts.
 */
template<
    typename Str1,
    typename Str2
>
unsigned
rev_comp_shifted_edit_distance( const Str1 & s1, const Str2 & s2 ) {
    BOOST_ASSERT( boost::size( s1 ) == boost::size( s2 ) );
    return std::min( shifted_edit_distance( s1, s2 ), shifted_edit_distance( make_rev_comp_range( s1 ), s2 ) );
}




/**
 * Calculate unique value for a sequence (shared with reverse complement).
 */
template< size_t alphabet_size, typename Seq >
seq_ordinal_t
calculate_seq_ordinal( const Seq & seq ) {
	return std::min(
		calculate_seq_ordinal_one_way< alphabet_size >( seq ),
		calculate_seq_ordinal_one_way< alphabet_size >( make_rev_comp_range( seq ) )
	);
}



/**
 * Make a geometric progression.
 */
template< typename OutputIt >
void
make_geometric_progression( int min_val, int max_val, double factor, OutputIt output_it ) {
	int last_val = 0;
	for( int curr_val = min_val; curr_val < max_val; curr_val = ( int )( curr_val * factor ) ) {
		// make sure we always increment by at least one.
		if( last_val == curr_val ) {
			if( max_val == ++curr_val ) {
				break;
			}
		}
		*output_it = curr_val;
		++output_it;
		// remember what we just added
		last_val = curr_val;
	}
	// Add the terminal value
	*output_it = max_val;
	++output_it;
}




/**
 * Find the best starting points.
 */
template<
	find_best_w_mers_storage_enum  Storage  = DEFAULT_FIND_WMER_STORAGE,
	typename                       Spec     = default_spec
>
struct start_finder
: boost::noncopyable
{
	STEME_TYPEDEF_SEQAN_TYPES( Spec )
	typedef data< Spec >                                data_t;                 ///< Data type.
	typedef model< Spec >                               model_t;                ///< Model type.
    typedef significance< Spec >                        significance_t;         ///< The significance type.
	typedef find_best_w_mers< Storage, Spec >          best_w_mers_finder;     ///< Type of find best W-mers.
	typedef typename best_w_mers_finder::ptr           best_w_mers_ptr;        ///< Pointer to best W-mers structure.
    typedef w_mer_evaluation< Spec >                    w_mer_eval_t;
    typedef typename w_mer_eval_t::vec                 eval_vec;
    typedef typename w_mer_eval_t::vec_ptr             eval_vec_ptr;

	/** A starting point for EM. */
	struct start
	: boost::less_than_comparable< start >
	, boost::equality_comparable< start >
	{
		double                                            score;                ///< The start's score.
		size_t                                             num_sites;            ///< The number of sites used in the start.
		std::string                                        w_mer;                ///< W-mer that seeds the start.
		eval_vec_ptr                                       best_w_mers;          ///< The best W-mers for the seed used to align the start.
		typename model_t::ptr                             _model;               ///< The model used to evaluate the start.

		/** Default constructor. */
		start()
		: score( 0.0 )
		, num_sites( 0 )
		{
		}

		/** Constructor. */
		start(
		    double score,
		    size_t num_sites,
		    std::string w_mer,
		    eval_vec_ptr best_w_mers,
		    typename model_t::ptr _model
		)
		: score( score )
		, num_sites( num_sites )
		, w_mer ( w_mer )
		, best_w_mers( best_w_mers )
		, _model( _model )
		{
		}

	    /** Comparison operator. */
	    bool
	    operator<( start other ) const {
	        if( this->score < other.score ) return true;
	        if( this->score > other.score ) return false;
	        if( this->num_sites < other.num_sites ) return true;
	        if( this->num_sites > other.num_sites ) return false;
	        if( this->w_mer < other.w_mer ) return true;
	        return false;
	    }

	    /// Equality comparison.
	    bool
	    operator==( start other ) const {
	        return this->score == other.score && this->num_sites == other.num_sites && this->w_mer == other.w_mer;
	    }
	};


	/**
	 * Descends the suffix tree to count the number of starts of a certain width. Ignores reverse complements in the
	 * count.
	 */
	struct start_counter_for_width {
		data_t &                           _data;               ///< The data.
		bool                               using_start_hash;    ///< True iff we are using a hash to remember which seeds (or reverse complements) we've already tried
		size_t                             W;                   ///< Width of starts we're counting.
		size_t                             count;               ///< The count of starts.
		seq_ordinal_seq                    starts_already_seen; ///< Remembers those starts we have seen already.

		/// Constructor.
		start_counter_for_width( data_t & _data, size_t W )
		: _data( _data )
		, using_start_hash( ordinal_type_large_enough< 4, seq_ordinal_t >( W ) )
		, W( W )
		{ }

		/// Count the starts in the index.
		size_t
		operator()() {
			count = 0;
			starts_already_seen.clear();
			visit_node( top_down_it( _data.index ) );
			return count;
		}

	protected:
		/// Visit a node in the index.
		void
		visit_node( top_down_it it ) {
			using namespace seqan;

			// ignore if N in parent edge
			if( _data.contains_n( W, it ) ) {
				return;
			}

			const size_t w = std::min( repLength( it ), W );

			// did we get to the bottom?
			if( W == w ) {

				// yes, we are at bottom

				// have we seen this W-mer (or its reverse complement) before?
				bool not_seen_w_mer_already = true;
				const prefix_t pre = prefix( representative( it ), w );
				if( using_start_hash ) {
					// if using hash, then actually check
					const seq_ordinal_t seq_ordinal = calculate_seq_ordinal< 4 >( pre );
					not_seen_w_mer_already = starts_already_seen.insert( seq_ordinal ).second;
				}
				if( not_seen_w_mer_already ) {
					++count;
				}

			} else {

				// continue : so go down and across
				if( goDown( it ) ) {
					visit_node( it );
					while( goRight( it ) ) {
						visit_node( it );
					}
				}
			}
		}
	};

	typedef boost::shared_ptr< start_finder >                    ptr;                   ///< Shared pointer.
	typedef size_t                                               partition_key;         ///< Key to index a partition.
	typedef std::vector< start >                                 start_vec;             ///< Vector of starts.
	typedef boost::shared_ptr< start_vec >                       start_vec_ptr;         ///< Pointer to vector of starts.
	typedef std::map< size_t, start_vec_ptr >                    start_map;             ///< Map from num_sites into vectors of best starts.
	typedef typename start_map::iterator                         start_map_it;          ///< Map from num_sites into best starts.
	typedef boost::shared_ptr< start_map >                       start_map_ptr;         ///< Pointer to partitions of starts.
	typedef std::vector< size_t >                                num_sites_vec;         ///< Vector of candidate num_sites.
	typedef boost::function<
		void (
		    eval_vec_ptr best_w_mers,
			const std::string & start,
			int num_sites,
			double score
		)
	>                                                        callback;              ///< Callback function.
	typedef std::vector< callback >                          callback_vec;          ///< Vector of callbacks.

	data_t &                           _data;               ///< Data.
	model_t &                          _model;              ///< Model.
    significance_t &                   _significance;       ///< The significance.
	size_t                             W;                   ///< Width of the model.
	num_sites_vec                      candidate_num_sites; ///< Candidates to be used as the number of sites.
	size_t                             num_to_find;         ///< Number of starts to find.
	start_map_ptr                      best_starts;         ///< The best starts found.
	seq_ordinal_seq                    starts_already_seen; ///< Remembers those starts we have seen already.
	bool                               using_start_hash;    ///< True iff we are using a hash to remember which seeds (or reverse complements) we've already tried
	efficiency_stats                   stats;               ///< Efficiency statistics for the descent of the suffix tree.
	size_t                             start_counter;       ///< Count the starts so we can ignore some of them.
	size_t                             starts_examined;     ///< How many starts we have examined.
	size_t                             speed_up;            ///< If >0, we only examine every n'th start.
	size_t                             first_to_examine;    ///< The index of the first start to examine.
	callback_vec                       cbs;                 ///< Callbacks for when a start is examined.


	/** Constructor. */
	start_finder(
		data_t &               _data,
		model_t &              _model,
        significance_t &       _significance,
		size_t                 min_num_sites,
		size_t                 max_num_sites,
		double                 candidate_starts_factor,
        size_t                 num_to_find = 1,
		size_t                 speed_up = 0,
		size_t                 first_to_examine = 0
	)
	: _data( _data )
	, _model( _model )
    , _significance( _significance )
	, W( _model.W() )
	, num_to_find( num_to_find )
	, best_starts( new start_map )
	, using_start_hash( ordinal_type_large_enough< 4, seq_ordinal_t >( _model.W() ) )
	, stats( _model.W() )
	, start_counter( 0 )
	, starts_examined( 0 )
	, speed_up( speed_up )
	, first_to_examine( first_to_examine )
	{
		if( speed_up && first_to_examine >= speed_up ) {
			throw std::logic_error( "Nonsensical arguments: first_to_examine >= speed_up");
		}
		make_geometric_progression(
			min_num_sites,
			max_num_sites,
			candidate_starts_factor,
			std::back_insert_iterator< num_sites_vec >( candidate_num_sites )
		);
	}


	/**
	 * Add a callback for when a start is examined.
	 */
	void
	register_callback( callback cb ) {
		cbs.push_back( cb );
	}


	/**
	 * Descend tree and find best starting points.
	 */
	void
	find_starts() {
	    // initialise statistics and counters
		stats.reset();
		start_counter = starts_examined = 0;

		// find the starts
		visit_node( top_down_it( _data.index ) );

		// sort the starts so that the best starts are first
		BOOST_FOREACH( start_vec_ptr starts, *best_starts | boost::adaptors::map_values ) {
		    using boost::lambda::bind;
		    boost::sort( *starts, bind( &start::score, boost::lambda::_1 ) > bind( &start::score, boost::lambda::_2 ) );
		}
	}


	/**
	 * Find starts just for given seed.
	 */
	template< typename Prefix >
	void
	find_starts_for_seed( const Prefix & pre ) {
		stats.reset();
		examine_start( pre );
	}



	/**
	 * Implementation to descend tree.
	 */
	void
	visit_node( top_down_it it ) {
		using namespace seqan;

		MS_DEBUG_STRING_VAR( w_mer, representative( it ) );

		// ignore if N in parent edge
		if( _data.contains_n( W, it ) ) {
			return;
		}

		const size_t w = std::min( repLength( it ), W );

		// did we get to the bottom?
		if( W == w ) {

			// yes, we are at bottom

			// have we seen this W-mer (or its reverse complement) before?
			bool not_seen_w_mer_already = true;
			const prefix_t pre = prefix( representative( it ), w );
			if( using_start_hash ) {
				// if using hash, then actually check
				const seq_ordinal_t seq_ordinal = calculate_seq_ordinal< 4 >( pre );
				not_seen_w_mer_already = starts_already_seen.insert( seq_ordinal ).second;
			}
			if( not_seen_w_mer_already ) {
				examine_start( pre );
			}

		} else {

			// continue : so go down and across
			if( goDown( it ) ) {
				visit_node( it );
				while( goRight( it ) ) {
					visit_node( it );
				}
			}
		}
	}



	eval_vec_ptr
	find_best_w_mers_for_model() {
        // find the best W-mers under the newly seeded model
        best_w_mers_ptr best_w_mer_finder = best_w_mers_ptr(
            new best_w_mers_finder( _data, _model, *candidate_num_sites.rbegin() )
        );
        best_w_mer_finder->descend_tree();
        stats += best_w_mer_finder->stats; // update start finder stats with best_w_mer_finder stats
        return best_w_mer_finder->best_w_mers;
	}


	eval_vec_ptr
    find_best_instances_for_model() {
        using boost::lambda::bind;

        // find the best instances under the newly seeded model
        find_instances< Spec > instance_finder( _data, _model, .2 );
        instance_finder.descend_tree();
        // sort them in decreasing Z order
        boost::sort( *instance_finder.instances, bind( &w_mer_eval_t::Z, boost::lambda::_1 ) > bind( &w_mer_eval_t::Z, boost::lambda::_2 ) );
        stats += instance_finder.stats; // update start finder stats with stats
        return instance_finder.instances;
    }


    inline
    bool
    is_close_to( const std::string & pre1, const std::string & pre2 ) {
        return rev_comp_shifted_edit_distance( pre1, pre2 ) <= boost::size( pre1 ) / 3;
    }


    /**
     * Do the instances overlap?
     *
     * \todo Make the threshold a parameter of the start finder.
     */
    inline
    bool
    do_instances_overlap( const eval_vec & instances1, const eval_vec & instances2 ) {
        const unsigned overlap = calculate_overlap( instances1.begin(), instances1.end(), W, instances2.begin(), instances2.end(), W );
        const double threshold = .3; // if about a third of the bases overlap we say the instances overlap
        return instances1.size() * W * threshold <= overlap;
    }


	template< typename Prefix >
	void
	examine_start( const Prefix & pre ) {
		using namespace seqan;

		MYRRH_ASSERT( size_t( boost::size( pre ) ) == W );
		const std::string prefix( std_string_from_seqan( pre ) );

		/// if we are ignoring some starts, then check if we should ignore this one.
		const bool ignore_this_start = speed_up && ( start_counter % speed_up != first_to_examine );
		++start_counter; //increment counter
		if( ! ignore_this_start ) {

			// update our count
			++starts_examined;

			// seed the model with the W-mer using pseudo-counts
			_model.bs.seed( pre, true );

			// find the best W-mers
			eval_vec_ptr instances = find_best_w_mers_for_model();
			BOOST_ASSERT( are_instances_sorted_by_score( *instances ) );

            // if we got enough sites
            if( instances->size() >= candidate_num_sites[ 0 ] ) {

                // for each candidate number of sites
                BOOST_FOREACH( size_t num_sites, candidate_num_sites ) {

                    // update the model using the best W-mers
                    const size_t actual_num_sites = _model.update_with_W_mers( *instances, num_sites );

                    // calculate the score for the model
                    const double score = - _significance.log_E_value( _model );

                    // get the vector of starts for this number of sites
                    start_vec_ptr & starts_ptr = ( *best_starts )[ num_sites ];
                    if( ! starts_ptr ) {
                        starts_ptr.reset( new start_vec );
                    }
                    start_vec & starts = *starts_ptr;

                    // do we have enough in the vector?
                    if( starts.size() < num_to_find ) {
                        // no, so add this start
                        eval_vec_ptr start_instances(
                            new eval_vec( instances->begin(), instances->begin() + actual_num_sites )
                        );
                        sort_instances_by_position( *start_instances );
                        starts.push_back(
                            start(
                                score,
                                actual_num_sites,
                                prefix,
                                start_instances,
                                _model.copy()
                            )
                        );
                    } else {
                        //find the worst start and see if this start is better
                        using boost::lambda::bind;
                        typename start_vec::iterator worst = boost::min_element( starts, bind( &start::score, boost::lambda::_1 ) < bind( &start::score, boost::lambda::_2 ) );
                        if( score > worst->score ) {
                            // This start is better than the worst saved so far.
                            //
                            // Now we need a vector of just the instances for this start
                            // so we can check the overlaps with other starts
                            // the instances need to be sorted by position as well
                            eval_vec_ptr start_instances(
                                new eval_vec( instances->begin(), instances->begin() + actual_num_sites )
                            );
                            sort_instances_by_position( *start_instances );


                            // We will replace the worst unless we can find
                            // another start that is close to this one.
                            //
                            // If there is another start close to this one,
                            // we will replace it unless it is better than
                            // this start when we will drop this start
                            //
                            typename start_vec::iterator i = starts.begin();
                            for( ; starts.end() != i; ++i ) {
                                if( worst != i && do_instances_overlap( *start_instances, *i->best_w_mers ) ) {
                                    if( score < i->score ) {
                                        // i is a close start that is better than this one
                                        // we are going to drop this start
                                        break;
                                    } else {
                                        // i is a close start that is worse than this one
                                        // is it worse than our current worst?
                                        if( worst->score < i->score ) {
                                            worst = i;
                                        }
                                    }
                                }
                            }
                            if( starts.end() == i ) { // if we did not break out of loop above
                                // update the start
                                *worst = start(
                                    score,
                                    actual_num_sites,
                                    prefix,
                                    start_instances,
                                    _model.copy()
                                );
                            }
                        }
                    }

                    // Make any callbacks we have.
                    BOOST_FOREACH( callback cb, cbs ) {
                        cb( instances, prefix, actual_num_sites, score );
                    }

                    // if we didn't have enough sites then don't bother trying for more.
                    if( actual_num_sites != num_sites ) {
                        // we did not have enough sites.
                        break;
                    }
                }
            }
		}

	}
};




} // namespace steme

#endif /* STEME_JR_13AUG2011_START_FINDER_H_ */
