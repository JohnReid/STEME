/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the type that finds the best W-mers in the STEME algorithm.
 *
 */

#ifndef STEME_JR_13AUG2011_FIND_BEST_W_MERS_H_
#define STEME_JR_13AUG2011_FIND_BEST_W_MERS_H_

#include <steme/descender.h>
#include <steme/descender_preferred_order.h>
#include <steme/w_mer_evaluation.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include <myrrh/store_n_best.h>

namespace steme {
/// Enumeration of the storage policys.
enum find_best_w_mers_storage_enum {
    store_using_multi_index,      ///< Store using a container indexed by score and position.
    store_in_set,                 ///< Store in a set.
    store_in_sorted_vector,       ///< Store in a sorted vector.
    store_in_sorted_list,         ///< Store in a sorted list.
    store_in_heap,                ///< Store in a heap.
    store_all                     ///< Store all.
};
#define DEFAULT_FIND_WMER_STORAGE ::steme::store_in_heap ///< Default storage type.



/**
 * Base class for best W-mers storage.
 */
template<
    typename Derived,
    typename Container,
    typename Spec = default_spec
>
struct find_best_w_mers_storage_base {
    typedef Container                                            container;         ///< Container.
    typedef typename container::iterator                         iterator;          ///< Iterator type.
    typedef typename container::const_iterator                   const_iterator;    ///< Const iterator type.
    typedef typename steme_seqan_types< Spec >::global_pos_t    global_pos_t;      ///< Index a global position.

    container &         c;                ///< The storage for the best W-mers.
    size_t              max_to_store;     ///< Maximum to store.
    size_t              num_stored;       ///< Maintain a count of the W-mers stored so far.

    /// Constructor.
    find_best_w_mers_storage_base( container & c, size_t max_to_store )
    : c( c )
    , max_to_store( max_to_store )
    { }

    /// Have we max'ed out the storage
    inline
    bool
    stored_enough() const {
        return num_stored == max_to_store;
    }

    /// The derived object.
    inline
    Derived &
    derived() {
        return *static_cast< Derived * >( this );
    }

    /// Reset the storage.
    inline
    void
    initialise() {
        c.clear();
        num_stored = 0;
        derived().reserve( max_to_store );
    }


    /// Accept one W-mer, returns if was inserted
    inline
    bool
    operator()( size_t W, global_pos_t global_pos, bool rev_comp, double Z ) {
        // work out where it would be inserted
        // \todo Don't need to construct the eval here, just check where
        // an evaluation with the given Z would be inserted
        w_mer_evaluation< Spec > eval( Z, global_pos, rev_comp );
        iterator insertion_point = derived().insertion_point( eval );

        // is there a better W-mer that overlaps?
        if( ! derived().does_better_overlap( W, global_pos, Z, insertion_point ) ) {
            // no - so insert the W-mer
            derived().insert( insertion_point, eval );
            return true;
        }

        return false;
    }
};



/**
 * Best W-mers storage sorted by Z mix-in.
 */
template<
    typename Derived,
    typename Container,
    typename Spec = default_spec
>
struct find_best_w_mers_storage_Z_sorted
{
    typedef typename Container::iterator         it;          ///< Iterator type.
    typedef typename steme_seqan_types< Spec >::global_pos_t    global_pos_t;      ///< Index a global position.

    /// Constructor.
    find_best_w_mers_storage_Z_sorted()
    { }

    /// for every better W-mer, check if it overlaps
    inline
    bool
    does_better_overlap( size_t W, global_pos_t global_pos, double Z, it wmer ) const {
        BOOST_FOREACH(
            const w_mer_evaluation< Spec > & better,
            static_cast< const Derived * >( this )->better_wmers( wmer )
        ) {
            if( wmer_overlap( better.global_pos, global_pos, W ) ) {
                return true;
            }
        }
        return false;
    }

    /// Remove any overlapping W-mers and discard worst W-mers.
    /// \return Worst Z
    inline
    double
    finalise( size_t W, size_t new_size ) {
#ifdef STEME_LOG_TIMING_INFO
        std::cout << "Got " << static_cast< const Derived * >( this )->c.size() << " W-mers at end of descent.\n";
#endif // STEME_LOG_TIMING_INFO

        if( ! static_cast< const Derived * >( this )->c.empty() ) {
            it i = static_cast< const Derived * >( this )->c.begin();
            ++i;
            size_t count = 1;
            while( count < new_size && static_cast< const Derived * >( this )->c.end() != i ) {
                if( does_better_overlap( W, i->global_pos, i->Z, i ) ) {
                    // this W-mer overlaps a better one so erase it.
                    i = static_cast< Derived * >( this )->erase( i );
                } else {
                    // this W-mer does not overlap a better one so keep it.
                    ++i;
                    ++count;
                }
            }
            static_cast< const Derived * >( this )->c.erase( i, static_cast< const Derived * >( this )->c.end() );
        }

        return static_cast< Derived * >( this )->get_worst();
    }
};



/**
 * Storage that holds best W-mers in a vector sorted by Z.
 */
template< typename Spec >
struct find_best_w_mers_sorted_vector_storage
: find_best_w_mers_storage_base<
    find_best_w_mers_sorted_vector_storage< Spec >,
    std::vector< w_mer_evaluation< Spec > >,
    Spec
>,
find_best_w_mers_storage_Z_sorted<
    find_best_w_mers_sorted_vector_storage< Spec >,
    std::vector< w_mer_evaluation< Spec > >,
    Spec
>
{
    /// Base type.
    typedef find_best_w_mers_storage_base<
        find_best_w_mers_sorted_vector_storage< Spec >,
        std::vector< w_mer_evaluation< Spec > >,
        Spec
    > base_t;

    typedef typename base_t::container container;
    typedef typename base_t::iterator iterator;

    /// Constructor.
    find_best_w_mers_sorted_vector_storage( container & c, size_t max_to_store )
    : base_t( c, max_to_store )
    { }

    /// Reserve a certain amount of space.
    inline
    void
    reserve( size_t space ) {
        base_t::c.reserve( space );
    }

    /// Get a bound on where to insert the evaluation.
    inline
    iterator
    insertion_point( w_mer_evaluation< Spec > eval ) {
        return std::upper_bound( base_t::c.begin(), base_t::c.end(), eval, std::greater< w_mer_evaluation< Spec > >() );
    }

    /// Get a range containing those W-mers that are better than the given W-mer
    inline
    boost::iterator_range< iterator >
    better_wmers( iterator wmer ) const {
        return boost::make_iterator_range( base_t::c.begin(), wmer );
    }

    /// Insert the evaluation
    /// \return Worst Z
    inline
    double
    insert( iterator insertion_point, const w_mer_evaluation< Spec > & eval ) {
        if( base_t::num_stored == base_t::max_to_store ) {
            base_t::c.resize( base_t::max_to_store - 1 );
            base_t::c.insert( insertion_point, eval );
        } else {
            base_t::c.insert( insertion_point, eval );
            ++base_t::num_stored;
        }
        return get_worst();
    }

    /// Erase the evaluation
    /// \return An iterator pointing at the W-mer after the erased one.
    inline
    iterator
    erase( iterator wmer ) {
        base_t::c.erase( wmer );
        return wmer;
    }

    /// Get the worst value in the storage
    inline
    double
    get_worst() const {
        return base_t::c.rbegin()->Z;
    }
};


/**
 * Storage that holds best W-mers in a list sorted by Z.
 */
template< typename Spec >
struct find_best_w_mers_sorted_list_storage
: find_best_w_mers_storage_base<
    find_best_w_mers_sorted_list_storage< Spec >,
    std::list< w_mer_evaluation< Spec > >,
    Spec
>,
find_best_w_mers_storage_Z_sorted<
    find_best_w_mers_sorted_list_storage< Spec >,
    std::list< w_mer_evaluation< Spec > >,
    Spec
>
{
    /// Base type.
    typedef find_best_w_mers_storage_base<
        find_best_w_mers_sorted_list_storage< Spec >,
        std::list< w_mer_evaluation< Spec > >
    > base_t;

    typedef typename base_t::container container;
    typedef typename base_t::iterator iterator;

    /// Constructor.
    find_best_w_mers_sorted_list_storage( container & c, size_t max_to_store )
    : base_t( c, max_to_store )
    { }

    /// Reserve a certain amount of space.
    inline
    void
    reserve( size_t space ) {
    }

    /// Get a bound on where to insert the evaluation.
    inline
    iterator
    insertion_point( w_mer_evaluation< Spec > eval ) {
        return std::upper_bound( base_t::c.begin(), base_t::c.end(), eval, std::greater< w_mer_evaluation< Spec > >() );
    }

    /// Get a range containing those W-mers that are better than the given W-mer
    inline
    boost::iterator_range< iterator >
    better_wmers( iterator wmer ) const {
        return boost::make_iterator_range( base_t::c.begin(), wmer );
    }

    /// Insert the evaluation
    /// \return Worst Z
    inline
    double
    insert( iterator insertion_point, const w_mer_evaluation< Spec > & eval ) {
        if( base_t::num_stored == base_t::max_to_store ) {
            base_t::c.resize( base_t::max_to_store - 1 );
            base_t::c.insert( insertion_point, eval );
        } else {
            base_t::c.insert( insertion_point, eval );
            ++base_t::num_stored;
        }
        return get_worst();
    }

    /// Erase the evaluation
    /// \return An iterator pointing at the W-mer after the erased one.
    inline
    iterator
    erase( iterator wmer ) {
        base_t::c.erase( wmer++ );
        return wmer;
    }

    /// Get the worst value in the storage
    inline
    double
    get_worst() const {
        return base_t::c.rbegin()->Z;
    }
};


/**
 * Storage that holds best W-mers in a set sorted by Z.
 */
template< typename Spec >
struct find_best_w_mers_set_storage
: find_best_w_mers_storage_base<
    find_best_w_mers_set_storage< Spec >,
    std::set< w_mer_evaluation< Spec >, std::greater< w_mer_evaluation< Spec > > >,
    Spec
>,
find_best_w_mers_storage_Z_sorted<
    find_best_w_mers_set_storage< Spec >,
    std::set< w_mer_evaluation< Spec >, std::greater< w_mer_evaluation< Spec > > >,
    Spec
>
{
	/// Base type.
	typedef find_best_w_mers_storage_base<
		find_best_w_mers_set_storage< Spec >,
		std::set< w_mer_evaluation< Spec >, std::greater< w_mer_evaluation< Spec > > >,
		Spec
	> base_t;

    typedef typename base_t::container container;
    typedef typename base_t::iterator iterator;

	/// Constructor.
	find_best_w_mers_set_storage( container & c, size_t max_to_store )
	: base_t( c, max_to_store )
	{ }

	/// Reserve a certain amount of space.
	inline
	void
	reserve( size_t space ) {
	}

	/// Get a bound on where to insert the evaluation.
	inline
	iterator
	insertion_point( w_mer_evaluation< Spec > eval ) {
		return base_t::c.upper_bound( eval );
	}

	/// Get a range containing those W-mers that are better than the given W-mer
	inline
	boost::iterator_range< iterator >
	better_wmers( iterator wmer ) const {
		return boost::make_iterator_range( base_t::c.begin(), wmer );
	}

	/// Insert the evaluation
	/// \return Worst Z
	inline
	double
	insert( iterator insertion_point, const w_mer_evaluation< Spec > & eval ) {
		if( base_t::num_stored == base_t::max_to_store ) {
		    base_t::c.insert( insertion_point, eval );
			// remove worst
			iterator to_erase = base_t::c.end();
			--to_erase;
			base_t::c.erase( to_erase );
		} else {
		    base_t::c.insert( insertion_point, eval );
			++base_t::num_stored;
		}
		return get_worst();
	}

	/// Erase the evaluation
	/// \return An iterator pointing at the W-mer after the erased one.
	inline
	iterator
	erase( iterator wmer ) {
	    base_t::c.erase( wmer++ );
		return wmer;
	}

	/// Get the worst value in the storage
	inline
	double
	get_worst() const {
		return base_t::c.rbegin()->Z;
	}
};


/**
 * Storage that holds best W-mers sorted by position and score.
 */
template< typename Spec >
struct find_best_w_mers_store_all
{
public:
    /// Container type that is part of the interface.
    typedef std::vector< w_mer_evaluation< Spec > > container;
    typedef typename steme_seqan_types< Spec >::global_pos_t    global_pos_t;      ///< Index a global position.

    container &            c;                ///< The container we use to store the best W-mers as part of the interface.
    int                    max_to_store;     ///< Maximum number to store.
    double                 worst;            ///< The worst value in the storage.

public:
    /// Constructor.
    find_best_w_mers_store_all( container & c, int max_to_store )
    : c( c )
    , max_to_store( max_to_store )
    { }

    /// Have we max'ed out the storage
    inline
    bool
    stored_enough() const {
        return int( c.size() ) >= max_to_store;
    }

    /// Get the worst value in the storage
    inline
    double
    get_worst() {
        return worst;
    }

    /// Reset the storage.
    inline
    void
    initialise() {
        c.clear();
        worst = std::numeric_limits< double >::max();
    }

    /// Accept one W-mer, returns true iff it was inserted
    inline
    bool
    operator()( size_t W, global_pos_t global_pos, bool rev_comp, double Z ) {
        c.push_back( w_mer_evaluation< Spec >( Z, global_pos, rev_comp ) );
        if( worst > Z ) {
            worst = Z;
        }
        return true;
    }

    /// Remove any overlapping W-mers and discard worst W-mers.
    /// \return Worst Z
    inline
    double
    finalise( size_t W, size_t new_size ) {
#ifdef STEME_LOG_TIMING_INFO
        std::cout << "Got " << c.size() << " W-mers at end of descent.\n";
        boost::timer t;
#endif // STEME_LOG_TIMING_INFO
        boost::sort( c, typename w_mer_evaluation< Spec >::comp_score() );
#ifdef STEME_LOG_TIMING_INFO
        std::cout << "Sort took " << t.elapsed() << " seconds.\n";
        t.restart();
#endif // STEME_LOG_TIMING_INFO
        typename w_mer_evaluation< Spec >::vec _c;
        output_nonoverlapping_from_score_sorted( c.begin(), c.end(), W, new_size, std::back_inserter( _c ) );
        std::swap( c, _c );
#ifdef STEME_LOG_TIMING_INFO
        std::cout << "Erase overlapping took " << t.elapsed() << " seconds.\n";
#endif // STEME_LOG_TIMING_INFO
        //erase_overlapping_from_score_sorted( c, W, new_size );
        return c.empty() ? -std::numeric_limits< double >::max() : c.rbegin()->Z;;
    }
};



/**
 * Storage that holds best W-mers sorted by position and score.
 */
template< typename Spec >
struct find_best_w_mers_multi_index_storage
{
protected:
    struct score {}; ///< Tag for multi_index.
    struct position {}; ///< Tag for multi_index.

    typedef ::boost::multi_index::tag< score    > score_tag;    ///< Multi-index tag for score.
    typedef ::boost::multi_index::tag< position > position_tag; ///< Multi-index tag for position.
    typedef typename steme_seqan_types< Spec >::global_pos_t    global_pos_t;      ///< Index a global position.

    /// Compares using positions
    struct comp_position
    {
        // compare a position and an eval
        bool operator()( global_pos_t global_pos, const w_mer_evaluation< Spec > & eval ) const { return global_pos < eval.global_pos; }

        // compare an eval and a score
        bool operator()( const w_mer_evaluation< Spec > & eval, global_pos_t global_pos ) const { return eval.global_pos < global_pos; }
    };

    /// Compares using scores
    struct comp_score
    {
        // compare a score and an eval
        bool operator()( double Z, const w_mer_evaluation< Spec > & eval ) const { return Z < eval.Z; }

        // compare an eval and a score
        bool operator()( const w_mer_evaluation< Spec > & eval, double Z ) const { return eval.Z < Z; }
    };

    /// Container indexed by score and by position separately.
    typedef boost::multi_index_container<
        w_mer_evaluation< Spec >,
        ::boost::multi_index::indexed_by<
            ::boost::multi_index::ordered_non_unique<
                score_tag,
                ::boost::multi_index::member<
                    w_mer_evaluation< Spec >,
                    double,
                    &w_mer_evaluation< Spec >::Z
                >
            >, // w-mers by score
            ::boost::multi_index::ordered_non_unique<
                position_tag,
                ::boost::multi_index::member<
                    w_mer_evaluation< Spec >,
                    global_pos_t,
                    &w_mer_evaluation< Spec >::global_pos
                >
            >  // w-mers by position
        >
    > mi_container;

    typedef typename mi_container::template index< score    >::type   wmers_by_score;    ///< W-mers indexed by score.
    typedef typename mi_container::template index< position >::type   wmers_by_position; ///< W-mers indexed by position.
    typedef typename wmers_by_position::iterator                       score_iterator;    ///< Iterator by score.
    typedef typename wmers_by_position::iterator                       position_iterator; ///< Iterator by position.

public:
    /// Container type that is part of the interface.
    typedef std::vector< w_mer_evaluation< Spec > > container;

protected:
    mi_container           mi;               ///< The multi-index container we use to store the best W-mers internally.
    container &            c;                ///< The multi-index container we use to store the best W-mers as part of the interface.
    int                    max_to_store;     ///< Maximum to store.
    int                    num_stored;       ///< Maintain a count of the W-mers stored so far.

public:
    /// Constructor.
    find_best_w_mers_multi_index_storage( container & c, int max_to_store )
    : c( c )
    , max_to_store( max_to_store )
    { }

    /// Have we max'ed out the storage
    inline
    bool
    stored_enough() const {
        return num_stored == max_to_store;
    }

    /// Get the worst value in the storage
    inline
    double
    get_worst() {
        return mi.empty() ? -std::numeric_limits< double >::max() : sorted_by_score().rbegin()->Z;
    }

    /// Reset the storage.
    inline
    void
    initialise() {
        mi.clear();
        num_stored = 0;
    }

    /// Accept one W-mer, returns true iff it was inserted
    inline
    bool
    operator()( size_t W, global_pos_t global_pos, bool rev_comp, double Z ) {
        position_iterator insertion_point = sorted_by_position().upper_bound( global_pos );

        if( ! does_better_overlap( W, global_pos, Z, insertion_point ) ) {
            sorted_by_position().insert( insertion_point, w_mer_evaluation< Spec >( Z, global_pos, rev_comp ) );
            return true;
        }

        return false;
    }

    /// Remove any overlapping W-mers and discard worst W-mers.
    /// \return Worst Z
    inline
    double
    finalise( size_t W, size_t new_size ) {
        c.clear();
        c.reserve( num_stored );
        output_non_overlapping_instances_from_sorted( sorted_by_position(), W, std::back_inserter( c ) );
        sort_instances_by_score( c );
        if( new_size < c.size() ) {
            c.resize( new_size );
        }
        return c.empty() ? -std::numeric_limits< double >::max() : c.rbegin()->Z;;
    }


protected:
    /// W-mers by score
    inline
    wmers_by_score &
    sorted_by_score() {
        return mi.template get< score >();
    }

    /// W-mers by position
    inline
    wmers_by_position &
    sorted_by_position() {
        return mi.template get< position >();
    }

    /// Does a better one overlap?
    inline
    bool
    does_better_overlap( int W, global_pos_t global_pos, double Z, position_iterator insert_point ) {
        // check those after
        position_iterator i = insert_point;
        typedef typename w_mer_evaluation< Spec >::global_pos_t global_pos_t;
        while( sorted_by_position().end() != i && i->global_pos < global_pos + global_pos_t( W ) ) {
            if( i->Z >= Z ) {
                return true;
            }
            ++i;
        }
        // check those before
        if( sorted_by_position().begin() != insert_point ) {
            i = insert_point;
            --i;
            while( global_pos < i->global_pos + global_pos_t( W ) ) {
                if( i->Z >= Z ) {
                    return true;
                }
                if( sorted_by_position().begin() == i ) {
                    break;
                }
                --i;
            }
        }
        return false;
    }

};


/**
 * Storage that holds best W-mers sorted by position and score.
 */
template< typename Spec >
struct find_best_w_mers_heap_storage
{
protected:
public:
    typedef std::vector< w_mer_evaluation< Spec > > container; ///< Container type that is part of the interface.
    typedef typename steme_seqan_types< Spec >::global_pos_t    global_pos_t;      ///< Index a global position.

protected:
    typedef myrrh::store_n_best< container, typename w_mer_evaluation< Spec >::comp_score > best_n_t; ///< Stores best W-mers.

    best_n_t               best_n;           ///< Stores best W-mers

public:
    /// Constructor.
    find_best_w_mers_heap_storage( container & c, int max_to_store )
    : best_n( max_to_store, c )
    { }

    /// Have we max'ed out the storage
    inline
    bool
    stored_enough() const {
        return best_n.at_capacity();
    }

public:
    /// Get the worst value in the storage
    inline
    double
    get_worst() {
        return best_n._storage.empty() ? -std::numeric_limits< double >::max() : best_n.get_worst().Z;
    }

    /// Reset the storage.
    inline
    void
    initialise() {
        best_n._storage.clear();
    }

    /// Accept one W-mer, returns true iff it was inserted
    inline
    bool
    operator()( size_t W, global_pos_t global_pos, bool rev_comp, double Z ) {
        return best_n( w_mer_evaluation< Spec >( Z, global_pos, rev_comp ) );
    }

    /// Remove any overlapping W-mers and discard worst W-mers.
    /// \return Worst Z
    inline
    double
    finalise( size_t W, size_t new_size ) {
#ifdef STEME_LOG_TIMING_INFO
        boost::timer t;
#endif // STEME_LOG_TIMING_INFO
        boost::sort( best_n._storage, typename w_mer_evaluation< Spec >::comp_score() );
#ifdef STEME_LOG_TIMING_INFO
        std::cout << "Sort took " << t.elapsed() << " seconds.\n";
        t.restart();
#endif // STEME_LOG_TIMING_INFO
        erase_overlapping_from_score_sorted< Spec >( best_n._storage, W, new_size );
#ifdef STEME_LOG_TIMING_INFO
        std::cout << "Removing overlaps took " << t.elapsed() << " seconds.\n";
#endif // STEME_LOG_TIMING_INFO
        return get_worst();
    }
};


/// Specialised to choose overlap policy type.
template< find_best_w_mers_storage_enum Storage, typename Spec >
struct find_best_w_mers_storage_selector {
};


/// Specialisation to choose store all for find best W-mers.
template< typename Spec >
struct find_best_w_mers_storage_selector< store_all, Spec > {
    typedef find_best_w_mers_store_all< Spec > type; ///< Type that holds best W-mers.
};


/// Specialisation to choose set storage for find best W-mers.
template< typename Spec >
struct find_best_w_mers_storage_selector< store_in_set, Spec > {
    typedef find_best_w_mers_set_storage< Spec > type; ///< Type that holds best W-mers.
};


/// Specialisation to choose sorted vector storage for find best W-mers.
template< typename Spec >
struct find_best_w_mers_storage_selector< store_in_sorted_vector, Spec > {
    typedef find_best_w_mers_sorted_vector_storage< Spec > type; ///< Type that holds best W-mers.
};


/// Specialisation to choose sorted list storage for find best W-mers.
template< typename Spec >
struct find_best_w_mers_storage_selector< store_in_sorted_list, Spec > {
    typedef find_best_w_mers_sorted_list_storage< Spec > type; ///< Type that holds best W-mers.
};


/// Specialisation to choose heap storage for find best W-mers.
template< typename Spec >
struct find_best_w_mers_storage_selector< store_in_heap, Spec > {
    typedef find_best_w_mers_heap_storage< Spec > type; ///< Type that holds best W-mers.
};


/// Specialisation to choose multi-index storage for find best W-mers.
template< typename Spec >
struct find_best_w_mers_storage_selector< store_using_multi_index, Spec > {
    typedef find_best_w_mers_multi_index_storage< Spec > type; ///< Type that holds best W-mers.
};


///// Specialisation to choose heap storage for find best W-mers.
//template<>
//struct find_best_w_mers_storage_selector< store_in_heap > {
//	typedef find_best_w_mers_heap_storage type; ///< Type that holds best W-mers.
//};





/**
 * Finds the best w-mers given the model.
 *
 * Note that MEME has slightly different functionality. MEME finds the best W-mer in each sequence subject
 * to the non-overlapping constraint. Our approach just finds the best W-mers regardless of how many occur
 * in each sequence.
 *
 * \todo Update to use PSPs.
 *
 * \todo Update to use heap for storage and compare runtimes
 */
template<
	find_best_w_mers_storage_enum    Storage = DEFAULT_FIND_WMER_STORAGE,
	typename                        Spec = default_spec
>
struct find_best_w_mers
: descender_preferred_order< find_best_w_mers< Storage, Spec >, Spec >
{
	STEME_TYPEDEF_SEQAN_TYPES( Spec )
	typedef find_best_w_mers< Storage, Spec >                                    this_t;            ///< This type.
	typedef boost::shared_ptr< this_t >                                           ptr;               ///< Pointer.
    typedef tree_descender< this_t, Spec >                                        base_t;            ///< Base type.
	typedef typename find_best_w_mers_storage_selector< Storage, Spec >::type   storage_policy;    ///< Storage policy type.
	typedef typename storage_policy::container                                   container_t;       ///< Container type to store best W-mers.
	typedef boost::shared_ptr< container_t >                                      container_ptr;     ///< Pointer to best W-mers.
	typedef w_mer_evaluation< Spec >                                              wmer_eval_t;

	size_t                               num_to_find;      ///< Number to find.
	container_ptr                        best_w_mers;      ///< Best w_mers.
	storage_policy                       storage;          ///< Stores the best W-mers.
	double                              worst_Z;          ///< The worst Z so far.
	double                              llr_threshold;    ///< The log likelihood ratio threshold.


	/** Constructor. */
	find_best_w_mers(
		typename base_t::data_t &         _data,
		typename base_t::model_t &        _model,
		size_t                            num_to_find
	)
	: descender_preferred_order< find_best_w_mers< Storage, Spec >, Spec >( _data, _model )
	, num_to_find( num_to_find )
	, best_w_mers( new container_t )
	, storage( *best_w_mers, num_to_find * 5 / 2 )
	{
	}


	/**
	 * Called by tree_descender base class as part of CRTP.
	 */
	void
	on_descent_begin() {
		base_t::stats.reset();
		worst_Z = llr_threshold = -std::numeric_limits< double >::max();
		storage.initialise();
	}


	/**
	 * \return The threshold on log likelihood ratio.
	 */
	double
	get_llr_threshold() {
		if( storage.stored_enough() ) {
			return llr_threshold;
		} else {
			return -std::numeric_limits< double >::max();
		}
	}


	/**
	 * Update the worst Z and log Z threshold term.
	 */
	inline
	void
	update_worst_Z( double new_worst_Z ) {
		if( new_worst_Z != worst_Z ) {
			worst_Z = new_worst_Z;
			llr_threshold = std::log( new_worst_Z / ( 1. - new_worst_Z ) );
		}
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
        if( ! storage.stored_enough() || worst_Z < Z_strand ) {

            const bool inserted = storage( base_t::W, global_pos, RevComp, Z_strand );
            if( inserted ) {
                // update the worst Z
                update_worst_Z( storage.get_worst() );
            }

        }
	}


	/**
	 * Called by tree_descender base class as part of CRTP.
	 */
	void
	on_descent_end() {
#ifdef STEME_LOG_TIMING_INFO
        boost::timer t;
#endif // STEME_LOG_TIMING_INFO

	    // get rid of any overlaps
	    storage.finalise( base_t::W, num_to_find );

        // update the worst Z
        update_worst_Z( storage.get_worst() );

#ifdef STEME_LOG_TIMING_INFO
        std::cout << "find_best_w_mers::on_descent_end() took " << t.elapsed() << " seconds.\n";
#endif // STEME_LOG_TIMING_INFO
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


	/// Are there any overlapping W-mers?
	bool
	has_overlapping() const {
		// for each W-mer
		for( typename container_t::iterator i = best_w_mers->begin(); i != best_w_mers->end(); ++i ) {
		    // for each better one
	        for( typename container_t::iterator better = best_w_mers->begin(); i != better; ++better ) {
	            // does it overlap?
	            if( wmer_overlap( better->global_pos, i->global_pos, base_t::W ) ) {
	                return true;
	            }
	        }
		}
		return false;
	}


	/**
	 * Updates the model with the best so many found sites. Also updates lambda.
	 *
	 * @return: How many sites used to update the model, this can be less than
	 * requested if not enough were found.
	 */
	size_t
	update_model_with_W_mers( size_t num_sites, bool use_pseudo_counts = false ) {
	    return base_t::_model.update_with_W_mers( *best_w_mers, num_sites, use_pseudo_counts );
	}
};






} // namespace steme

#endif /* STEME_JR_13AUG2011_FIND_BEST_W_MERS_H_ */
