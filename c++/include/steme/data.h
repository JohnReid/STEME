/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the data type for STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_DATA_H_
#define STEME_JR_13AUG2011_DATA_H_

#include <steme/seqan_types.h>

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/shared_ptr.hpp>

namespace steme {



/// forward decl
template< typename Spec = default_spec >
struct data;






/**
 * An abstract base class for a base-pair resolution PSP.
 */
struct base_resolution_psp {
	typedef boost::shared_ptr< base_resolution_psp >    ptr;             ///< Pointer to position-specific prior.
	typedef long unsigned int                         pos_t;           ///< Type that identifies a position.

	/// Return the maximal PSP for all the bases in the given range
	virtual
	double
	get_maximal( pos_t seq, pos_t seq_pos, size_t W ) = 0;

	/// Adjust the PSP to account for a binding site at the given location with the given expected Z.
	virtual
	void
	adjust_for_binding_site( pos_t seq, pos_t seq_pos, size_t W, double Z ) = 0;
};


/**
 * Dense base-pair resolution priors.
 */
struct dense_base_resolution_psp
: public base_resolution_psp, boost::noncopyable
{
	typedef std::vector< double >                psp_vec;               ///< A sparse vector of doubles.
	typedef std::vector< psp_vec >               psp_vec_vec;           ///< A vector of sparse vectors of doubles.

	psp_vec_vec storage; ///< Stores the priors.

	/// Constructor.
	template< typename Index >
	dense_base_resolution_psp( Index & index ) {
		using namespace seqan;
		typedef typename seqan::Host< Index >::Type text_t;
		typedef typename seqan::Id< text_t >::Type id_t;

		text_t & text = getFibre( index, EsaText() );
		storage.resize( length( text ) );
		for( size_t i = 0; storage.size() != i; ++i ) {
			const id_t _id = positionToId( text, i );
			storage[i].resize( length( valueById( text, _id ) ) );
			std::fill( storage[ i ].begin(), storage[ i ].end(), 0. );
		}
	}

	virtual ~dense_base_resolution_psp() { }

	/// Return the maximal PSP for all the bases in the given range
	virtual
	double
	get_maximal( base_resolution_psp::pos_t seq, base_resolution_psp::pos_t seq_pos, size_t W ) {
		// find the highest position-specific prior in this W-mer. This will be the position least
		// likely to contain a binding site a priori
		const psp_vec & seq_priors = storage[ seq ];
		psp_vec::const_iterator W_mer_begin = seq_priors.begin() + seq_pos;
                const double & (*pt2Max)(const double &, const double &) = std::max< double >;
		const double g = std::accumulate(
			W_mer_begin,
			W_mer_begin + W,
			0.,
			pt2Max
		);
		return g;
	}

	/// Adjust the PSP to account for a binding site at the given location with the given expected Z.
	virtual
	void
	adjust_for_binding_site( base_resolution_psp::pos_t seq, base_resolution_psp::pos_t seq_pos, size_t W, double Z ) {
		psp_vec & seq_priors = storage[ seq ];
		psp_vec::iterator i = seq_priors.begin() + seq_pos;
		psp_vec::iterator W_mer_end = i + W;
		for( ; W_mer_end != i; ++i ) {
			*i = std::min( 1., *i + Z );
		}
	}
};


/**
 * Sparse base-pair resolution priors.
 */
struct sparse_base_resolution_psp
: public base_resolution_psp, boost::noncopyable
{
	typedef boost::numeric::ublas::mapped_vector< double > psp_vec;      ///< A sparse vector of doubles.
	typedef std::vector< psp_vec >                         psp_vec_vec;  ///< A vector of sparse vectors of doubles.

	psp_vec_vec storage; ///< Stores the priors.

	/// Constructor.
	template< typename Index >
	sparse_base_resolution_psp( Index & index ) {
		using namespace seqan;
		typedef typename seqan::Host< Index >::Type text_t;
		typedef typename seqan::Id< text_t >::Type id_t;

		text_t & text = getFibre( index, EsaText() );
		storage.resize( length( text ) );
		for( size_t i = 0; storage.size() != i; ++i ) {
			const id_t _id = positionToId( text, i );
			storage[i].resize( length( valueById( text, _id ) ) );
		}
	}

	virtual ~sparse_base_resolution_psp() { }

	/// Return the maximal PSP for all the bases in the given range
	virtual
	double
	get_maximal( base_resolution_psp::pos_t seq, base_resolution_psp::pos_t seq_pos, size_t W ) {
		// find the highest position-specific prior in this W-mer. This will be the position least
		// likely to contain a binding site a priori
		const psp_vec & seq_priors = storage[ seq ];
                const double & (*pt2Max)(const double &, const double &) = std::max< double >;
		const double g = std::accumulate(
			seq_priors.find( seq_pos ),
			seq_priors.find( seq_pos + W ),
			0.,
			pt2Max
		);
		return g;
	}

	/// Adjust the PSP to account for a binding site at the given location with the given expected Z.
	virtual
	void
	adjust_for_binding_site( base_resolution_psp::pos_t seq, base_resolution_psp::pos_t seq_pos, size_t W, double Z ) {
		psp_vec & seq_priors = storage[ seq ];
		for( size_t i = seq_pos; seq_pos + W != i; ++i ) {
			seq_priors[ i ] = std::min( 1., seq_priors[ i ] + Z );
		}
	}
};






/**
 * Abstract base class for position-specific priors for a model of a particular width.
 */
template< typename Spec >
struct psp {

    typedef boost::shared_ptr< psp >    ptr;             ///< Pointer to position-specific prior.
    typedef long unsigned int         pos_t;           ///< Type that identifies a position.
    STEME_TYPEDEF_SEQAN_TYPES( Spec )

    /// Return g at the given position.
    virtual
    double
    get_g( pos_t seq, pos_t seq_pos ) = 0;

    /// get the bound on g
    virtual
    double
    get_min_g( top_down_it it ) = 0;

    /// update after base-resolution priors have changed.
    virtual
    void
    update() = 0;
};




/**
 * Calculates model's priors from the base resolution priors in the data.
 *
 * Uses CRTP to store values in a derived class.
 */
template<
    typename Derived,
    typename Spec
>
struct model_priors_calculator {

    typedef data< Spec >                        data_t;             ///< The data type.

    data_t &             _data;              ///< The data.
    size_t               W;                  ///< The width of the motif we are calculating the priors for.
    STEME_TYPEDEF_SEQAN_TYPES( Spec )

    /// Constructor.
    model_priors_calculator( data_t & _data, size_t W )
    : _data( _data )
    , W( W )
    { }

    /// Calculate model's priors from the base resolution priors in the data.
    void
    calculate_weights() {
        visit_node( top_down_it( _data.index ) );
    }

protected:

    /// Visit a node to calculate its minimum and maximum
    double
    visit_node( top_down_it it ) {
        using namespace seqan;

        MS_DEBUG_STRING_VAR( w_mer, representative( it ) );

        // used to identify this node.
        typename VertexDescriptor< index_t >::Type v = value( it );

        // the min g in the sub-tree under this node
        double min_g;

        // if we haven't reached depth W in tree, visit the node's children and update the maximum weight from them
        const size_t w = repLength( it );
        if( w < W ) {
            if( goDown( it ) ) {
                min_g = visit_node( it );
                while( goRight( it ) ) {
                    const double child_min_g = visit_node( it );
                    if( child_min_g < min_g ) {
                        min_g = child_min_g;
                    }
                }
            } else { // we had no children and we're not deep enough in the tree so the minimum g is 1.
                min_g = 1.;
            }
        } else {
            // we are at depth W (or greater)
            // calculate the minimum g for any occurrence at this node
            min_g = 1.;
            occurrences_t occs = seqan::getOccurrences( it );
            typedef typename seqan::Size< occurrences_t >::Type occ_size_t;
            const occ_size_t num_occurrences = seqan::length( occs );
            for( occ_size_t i = 0; num_occurrences != i; ++i ) {

                // get this occurrence
                const typename Value< occurrences_t >::Type & occ = occs[ i ];

                // find the highest position-specific prior in this W-mer. This will be the position least
                // likely to contain a binding site a priori
                const double g = _data._base_prior->get_maximal( occ.i1, occ.i2, W );

                // store weight
                static_cast< Derived * >( this )->store_g( occ.i1, occ.i2, g );

                // update the bound
                if( g < min_g ) {
                    min_g = g;
                }
            }
        }

        // store minimum if it is not 0.
        if( min_g > 0. ) {
            // term will be infinite if min_g is 1.
            static_cast< Derived * >( this )->store_min_g( v, min_g );
        }

        // return the minimum g for this sub-tree
        return min_g;
    }
};




/**
 * Dense priors for a model of a particular width.
 */
template< typename Spec = default_spec >
struct dense_psp
: psp< Spec >
, boost::noncopyable
, model_priors_calculator< dense_psp< Spec >, Spec >
{
    typedef typename psp< Spec >::pos_t          pos_t;           ///< Index into a position.
    typedef data< Spec >                         data_t;          ///< The data type.
    typedef seqan::String< double >              vertex_g_vec;    ///< Vector of weights indexed by suffix array vertex.
    typedef std::vector< double >                pos_g_vec;       ///< Vector of weights
    typedef std::vector< pos_g_vec >             vec_pos_g_vecs;  ///< Vector of vector of weights.
    STEME_TYPEDEF_SEQAN_TYPES( Spec )

protected:
    vertex_g_vec                min_g;         ///< Minimum g for any X_n under a suffix tree node.
    vec_pos_g_vecs              g;             ///< The weights.

public:
    /** Constructor. Descend the suffix tree to calculate the weights from the position-specific priors. */
    dense_psp( data_t & _data, size_t W )
    : model_priors_calculator< dense_psp< Spec >, Spec >( _data, W )
    , g( _data.num_sequences() )
    {
        // resize our vectors of weights for each sequence
        for( size_t i = 0; g.size() != i; ++i ) {
            g[ i ].resize( _data.get_sequence_length( i ) );
        }

        // resize our map from vertices to minimum g
        seqan::resizeVertexMap( _data.index, min_g );

        update();
    }

    virtual ~dense_psp() { }

    /// update after base-resolution priors have changed.
    virtual
    void
    update() {
        // fill with 0s
        seqan::arrayFill( seqan::begin( min_g ), seqan::end( min_g ), 0. );

        // calculate the weights
        model_priors_calculator< dense_psp< Spec >, Spec >::calculate_weights();
    }

    /// Return the weight at the given position.
    virtual
    double
    get_g( pos_t seq, pos_t seq_pos ) {
        return g[ seq ][ seq_pos ];
    }

    /// Get a bound on the prior
    virtual
    double
    get_min_g( top_down_it it ) {
        return getProperty( min_g, value( it ) );
    }

    /// Store g
    inline
    void
    store_g( pos_t seq, pos_t seq_pos, double g ) {
        this->g[ seq ][ seq_pos ] = g;
    }

    /// Store the bound on g
    inline
    void
    store_min_g( const typename seqan::VertexDescriptor< index_t >::Type & v, double g ) {
        // assign the minimum prior for this sub-tree to the node
        assignProperty( min_g, v, g );
    }
};




/**
 * Sparse prior weights for a model of a particular width.
 */
template< typename Spec = default_spec >
struct sparse_psp
: psp< Spec >
, boost::noncopyable
, model_priors_calculator< sparse_psp< Spec >, Spec >
{
    typedef typename psp< Spec >::pos_t                       pos_t;           ///< Index into a position.
    typedef data< Spec >                                       data_t;          ///< The data type.
    typedef unsigned long                                     vertex_id_t;     ///< Type that identifies a vertex.
    typedef std::map< vertex_id_t, double >                    vertex_g_map;    ///< Map from vertices to weights.
    typedef std::map< pos_t, double >                          pos_g_map;       ///< Map from positions to weights.
    typedef std::vector< pos_g_map >                            vec_pos_g_maps;  ///< Vector of maps from positions to weights.
    STEME_TYPEDEF_SEQAN_TYPES( Spec )

protected:
    vertex_g_map         min_g;         ///< Minimum g for all X_n under this node.
    vec_pos_g_maps       g;             ///< The weights.

public:
    /** Constructor. Descend the suffix tree to calculate the weights from the position-specific priors. */
    sparse_psp( data_t & _data, size_t W )
    : model_priors_calculator< sparse_psp< Spec >, Spec >( _data, W )
    {
        update();
    }

    virtual ~sparse_psp() { }

    /// update after base-resolution priors have changed.
    virtual
    void
    update() {
        g.clear();
        g.resize( model_priors_calculator< sparse_psp< Spec >, Spec >::_data.num_sequences() );

        // calculate the weights
        model_priors_calculator< sparse_psp< Spec >, Spec >::calculate_weights();
    }

    /// Get a bound on the prior
    virtual
    double
    get_min_g( top_down_it it ) {
        using namespace seqan;

        // look for entry in our map
        vertex_g_map::iterator i = min_g.find( _getId( value( it ) ) );
        return min_g.end() == i ? 0. : i->second;
    }

    /// Return the weight at the given position.
    virtual
    double
    get_g( pos_t seq, pos_t seq_pos ) {
        pos_g_map & seq_map = g[ seq ];
        typename pos_g_map::iterator i = seq_map.find( seq_pos );
        return seq_map.end() == i
            ? 0.
            : i->second
            ;
    }

    /// Store g
    inline
    void
    store_g( pos_t seq, pos_t seq_pos, double g ) {
        // only store weights that are not 1.
        if( g != 0. ) {
            this->g[ seq ][ seq_pos ] = g;
        }
    }

    /// Store the bound on g
    inline
    void
    store_min_g( const typename seqan::VertexDescriptor< index_t >::Type & v, double bound ) {
        min_g[ seqan::_getId( v ) ] = bound;
    }
};


/**
 * Holds a concatenated string of all sequences to make access to W-mers quick.
 */
template< typename Derived, typename Spec >
struct wmer_concatenated_string_policy {
    STEME_TYPEDEF_SEQAN_TYPES( Spec )
    string_t        concatenated_seqs;  ///< A concatenation of all the strings in the index's text (this is a copy of the strings).

protected:
    void initialise() {
        // copy all the strings into one long string
        // this means get_W_mer() can be much more efficient.
        BOOST_FOREACH( string_t & seq, static_cast< Derived * >( this )->get_text() ) {
            seqan::append( concatenated_seqs, seq );
        }
    }

public:
    /**
     * Gets the W-mer at the given global position.
     *
     * Is faster than the previous method that used posLocalize to get to a local position first.
     */
    infix_t
    get_W_mer( size_t W, global_pos_t pos ) {
        return seqan::infixWithLength( concatenated_seqs, pos, W );
    }
};


/**
 * Accesses W-mers by converting a global position into a seq # and offset. Does not need to store a
 * concatenated string.
 */
template< typename Derived, typename Spec >
struct wmer_simple_policy {
    STEME_TYPEDEF_SEQAN_TYPES( Spec )

protected:
    /** No initialisation to do. */
    void initialise() { }

public:
    /**
     * Get the W-mer that starts at the given index.
     */
    infix_t
    get_W_mer( size_t W, global_pos_t pos ) {
        using namespace seqan;
        seq_index_t id;
        seq_offset_t offset;
        boost::tie( id, offset ) = static_cast< Derived * >( this )->pos_localise( pos );
#ifdef MS_DO_DEBUG_CHECKS
        if( offset < 0 || length( getValueById( static_cast< Derived * >( this )->get_text(), id ) ) < offset + W ) {
            throw std::logic_error(
                MS_MAKE_STRING(
                    "No W-mer at given position (" << pos << "). That sequence is too short. Id="
                    << id << ", offset=" << offset
                )
            );
        }
#endif //MS_DO_DEBUG_CHECKS
        return infix( getValueById( static_cast< Derived * >( this )->get_text(), id ), offset, offset+W );
    }
};



/**
 * Chooses the policy for how the data class handles retrieving W-mers.
 */
template< typename Spec >
struct choose_wmer_policy {
    typedef wmer_simple_policy< data< Spec >, Spec > type; ///< choose a concatenated string policy by default.
};



/**
 * Aggregates the data.
 *
 * \todo Implement sparse storage for PSP
 *
 * \todo Provide python interface to PSP
 */
template< typename Spec >
struct data
: boost::noncopyable
, choose_wmer_policy< Spec >::type
{
    typedef typename choose_wmer_policy< Spec >::type     wmer_policy_t;         ///< Mixin class that handles how to get W-mers.
	typedef data< Spec >                                   self_t;                ///< Self type.
	typedef boost::shared_ptr< self_t >                    ptr;                   ///< Pointer.
	typedef unsigned char                                first_unknown_t;       ///< Type to store position of first unknown.
	typedef seqan::String< first_unknown_t >               first_unknown_vec;     ///< Vector of sizes.
	typedef std::vector< size_t >                          _counts_t;             ///< Type to store counts for each width/depth.
	typedef std::map< size_t, typename psp< Spec >::ptr >  priors_cache_t;      ///< Type to cache priors for specific widths.
    STEME_TYPEDEF_SEQAN_TYPES( Spec )

	index_t &                            index;              ///< The index over our sequences.
	const raw_text_t &                   raw_text;           ///< A concatenation of all the strings in the index's text (this view is provided by the index).
	size_t                               N;                  ///< The number of bases in the sequences (i.e. number of global positions).
	first_unknown_vec                    first_unknown;      ///< Stores index of first position of unknown base in representative of each node.
    size_t                               max_W;              ///< Maximum width to consider for counts.
    _counts_t                            occurrence_counts;  ///< How many occurrences at each depth (excepting W-mers with unknown bases).
    bool                                 using_dense;        ///< Are we using dense priors or not?
	typename base_resolution_psp::ptr    _base_prior;        ///< Position-specific prior for each base. Prior probability base is not available for binding.
	priors_cache_t                       _priors_cache;      ///< Caches prior for each width.


	/** Constructor. */
	data( index_t & index, bool using_dense = false, size_t max_W = 0 )
		: index( index )
		, raw_text( seqan::getFibre( index, seqan::EsaRawText() ) )
		, N( seqan::length( index ) )
	    , max_W( max_W )
	    , occurrence_counts( max_W + 1 )
	    , using_dense( using_dense )
		, _base_prior(
			using_dense
				? base_resolution_psp::ptr( new  dense_base_resolution_psp( index ) )
				: base_resolution_psp::ptr( new sparse_base_resolution_psp( index ) )
		)
	{
	    if( max_W > std::numeric_limits< first_unknown_t >::max() ) {
	        throw std::invalid_argument( "max_W too big for first_unknown vector." );
	    }
	    wmer_policy_t::initialise();
		preprocess_index();
	}


	/// Descend the suffix tree and store the length of the representative of each node and whether it contains an 'N'.
	void
	preprocess_index() {
        seqan::resizeVertexMap( index, first_unknown );
		visit_node( top_down_it( index ), 0, std::numeric_limits< size_t >::max() );
	}


	/// Visit a node to calculate the length of its representative.
	void
	visit_node( top_down_it it, size_t last_length, size_t first_unknown_index ) {
		using namespace seqan;

		// get the new length.
        const size_t this_length = repLength( it );

        // index into our property vectors
        const typename Value< top_down_it >::Type idx = value( it );

        // look for an unknown base
        const int first_unknown_in_parent =
            find_first_unknown(
                infix(
                    representative( it ),
                    last_length,
                    max_W ? std::min( max_W, this_length ) : this_length
                )
            );
        if( -1 != first_unknown_in_parent ) {
            first_unknown_index = first_unknown_in_parent + last_length;
        }
        assignProperty( first_unknown, idx, first_unknown_index );

        // update the counts of occurrences at each depth until we reach the end of this representative, the maximum W or the first
        // unknown
        const size_t end = std::min( this_length + 1, first_unknown_index );
        if( 0 == max_W && end > occurrence_counts.size() ) {
            occurrence_counts.resize( end, 0 ); // resize if needed
        }
        const size_t num_known = countOccurrences( it );
        for( size_t w = last_length + 1; ( ! max_W || max_W >= w ) && end > w; ++w ) {
            // add the occurrence count
            BOOST_ASSERT( w < occurrence_counts.size() );
            occurrence_counts[ w ] += num_known;
        }

        // only go deeper if we haven't reached maximum depth already and we don't have an unknown base
        if( ( ! max_W || max_W > this_length ) && std::numeric_limits< size_t >::max() == first_unknown_index ) {
            // visit the node's children
            if( goDown( it ) ) {
                visit_node( it, this_length, first_unknown_index );
                while( goRight( it ) ) {
                    visit_node( it, this_length, first_unknown_index );
                }
            }
		}
	}


    /// Update the width-specific priors after the base-resolution prior have changed.
    void
    update_priors() {
        using boost::lambda::bind;
        using boost::adaptors::map_values;
        boost::for_each( _priors_cache | map_values, bind( &psp< Spec >::update, *boost::lambda::_1 ) );
    }

    /// Get the prior for the given width.
    typename psp< Spec >::ptr
    prior_for_width( size_t W ) {
        typename priors_cache_t::iterator i = _priors_cache.find( W );
        if( _priors_cache.end() == i ) {
            boost::tie( i, boost::tuples::ignore ) = _priors_cache.insert(
                typename priors_cache_t::value_type(
                    W,
                    using_dense
                        ? typename psp< Spec >::ptr( new  dense_psp< Spec >( *this, W ) )
                        : typename psp< Spec >::ptr( new sparse_psp< Spec >( *this, W ) )
                )
            );
        }
        return i->second;
    }


    /// Number of sequences in the data.
    size_t
    num_sequences() const {
        return seqan::length( seqan::getFibre( index, seqan::EsaText() ) );
    }


    /// Get the length of given sequence.
    size_t
    get_sequence_length( size_t sequence ) {
        using namespace seqan;
        text_t & text = get_text();
        const id_t _id = positionToId( text, sequence );
        return length( valueById( text, _id ) );
    }

    /** Number of W-mers in index. */
    size_t
    num_W_mers( size_t W ) {
        const size_t num_seqs = seqan::length( get_text() );
        return N - num_seqs * (W-1);
    }


    /// How many known occurrences at the given depth?
    size_t
    get_occurrence_count( size_t W ) {
        return W > occurrence_counts.size() ? 0 : occurrence_counts[ W ];
    }


	/// Is there an 'N' in the first W-mer of the representative?
    inline
	bool
	contains_n( size_t W, top_down_it it ) {
		return seqan::getProperty( first_unknown, value( it ) ) < W;
	}


	/// Return the length of the representative of the node.
	inline
	size_t
	get_rep_length( top_down_it it ) {
	    return seqan::repLength( it );
	}


	/// Get the string set limits that are used to convert between local and global positions.
	const string_set_limits_t &
	get_string_set_limits() {
		return seqan::stringSetLimits( seqan::getFibre( index, seqan::EsaText() ) );
	}


	/** Compare 2 iterators. */
	static
	bool
	compare_iterators( top_down_it i1, top_down_it i2 ) {
		return seqan::representative( i1 ) < seqan::representative( i2 );
	}


	/** An iterator at the root of the index. */
	top_down_it
	index_root() {
		return top_down_it( index );
	}


	/** Do 2 global positions overlap? */
	static
	bool
	overlap( size_t W , global_pos_t pos1, global_pos_t pos2 ) {
		return wmer_overlap( pos1, pos2, W );
	}


	/** Returns a tuple: sequence id, offset for a global position. */
	local_pos_t
	pos_localise( global_pos_t pos ) {
		using namespace seqan;
		sa_value_t local_pos;
		posLocalize( local_pos, pos, get_string_set_limits() );
		return boost::make_tuple(
			seqan::positionToId( get_text(), local_pos.i1 ),
			getSeqOffset( local_pos, get_string_set_limits() )
		);
	}


	/** Converts a local position to a global position. */
	global_pos_t
	pos_globalise( seq_index_t id, seq_offset_t offset ) {
		return seqan::posGlobalize( sa_value_t( id, offset ), get_string_set_limits() );
	}


	/// Get the index's text.
	text_t &
	get_text() {
		return seqan::getFibre( index, seqan::EsaText() );
	}
};

} // namespace steme


#endif /* STEME_JR_13AUG2011_DATA_H_ */
