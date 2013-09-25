/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the type that represents the evaluation of a W-mer under the model.
 *
 */

#ifndef STEME_JR_6SEP2011_W_MER_EVALUATION_H_
#define STEME_JR_6SEP2011_W_MER_EVALUATION_H_

#include <steme/defs.h>
#include <steme/seqan_types.h>

#include <boost/range/algorithm/sort.hpp>
#include <boost/icl/interval.hpp>


namespace steme {


/**
 * An evaluation of a W-mer.
 */
template< typename Spec = default_spec >
struct w_mer_evaluation
: boost::less_than_comparable< w_mer_evaluation< Spec > >
, boost::equality_comparable< w_mer_evaluation< Spec > >
{
    typedef std::vector< w_mer_evaluation< Spec > >     vec;        ///< A vector of W-mer evaluations.
    typedef boost::shared_ptr< vec >                     vec_ptr;    ///< A pointer to a vector of W-mer evaluations.
    STEME_TYPEDEF_SEQAN_TYPES( Spec )

    double          Z;                      ///< Z for this W-mer.
    global_pos_t     global_pos;             ///< Global position.
    bool            rev_comp;               ///< Was it the reverse complement that scored highly.

    /** Default constructor. */
    w_mer_evaluation()
    : Z( 0. )
    , global_pos( 0 )
    , rev_comp( false )
    { }

    /** Constructor. */
    w_mer_evaluation( double Z, global_pos_t global_pos, bool rev_comp )
    : Z( Z )
    , global_pos( global_pos )
    , rev_comp( rev_comp )
    { }

    /** Comparison operator. */
    bool operator<( w_mer_evaluation other ) const {
        if( this->Z < other.Z ) return true;
        if( this->Z > other.Z ) return false;
        if( this->global_pos < other.global_pos ) return true;
        if( this->global_pos > other.global_pos ) return false;
        if( this->rev_comp < other.rev_comp ) return true;
        return false;
    }

    /// Equality comparison.
    bool operator==( w_mer_evaluation other ) const {
        return this->Z == other.Z && this->global_pos == other.global_pos && this->rev_comp == other.rev_comp;
    }

    /// Compares using scores
    struct comp_score
    {
        // compare a score and an eval
        bool operator()( const w_mer_evaluation & lhs, const w_mer_evaluation & rhs ) const { return lhs.Z > rhs.Z; }
    };

    /// Compares using positions
    struct comp_position
    {
        // compare a score and an eval
        bool operator()( const w_mer_evaluation & lhs, const w_mer_evaluation & rhs ) const { return lhs.global_pos > rhs.global_pos; }
    };

};



/** Sort a range of instances by their position. */
template< typename Instances >
void
sort_instances_by_position( Instances & instances ) {
    using boost::lambda::bind;
    typedef typename boost::range_value< Instances >::type w_mer_eval_t;
    boost::sort( instances, bind( &w_mer_eval_t::global_pos, boost::lambda::_1 ) < bind( &w_mer_eval_t::global_pos, boost::lambda::_2 ) );
}


/** Sort a range of instances into descending scores. */
template< typename Instances >
void
sort_instances_by_score( Instances & instances ) {
    using boost::lambda::bind;
    typedef typename boost::range_value< Instances >::type w_mer_eval_t;
    boost::sort( instances, bind( &w_mer_eval_t::Z, boost::lambda::_1 ) > bind( &w_mer_eval_t::Z, boost::lambda::_2 ) );
}


/** Is the range of instances sorted into descending scores. */
template< typename Instances >
bool
are_instances_sorted_by_score( const Instances & instances ) {
    if( ! boost::empty( instances ) ) {
        typedef typename boost::range_iterator< const Instances >::type it;

        it i = boost::begin( instances );
        it j = boost::begin( instances );
        ++i;
        while( boost::end( instances ) != i ) {
            if( i->Z > j->Z ) {
                return false;
            }
            ++i;
            ++j;
        }
    }
    return true;
}


/** Output non-overlapping instances. */
template< typename Instances, typename OutIt >
void
output_non_overlapping_instances( Instances & instances, int W, OutIt out, bool already_sorted_by_position = false ) {
    if( ! already_sorted_by_position ) {
        sort_instances_by_position( instances );
    }
    output_non_overlapping_instances_from_sorted( instances, W, out );
}


/** Output non-overlapping instances given a list of instances sorted by position. */
template< typename Instances, typename OutIt >
void
output_non_overlapping_instances_from_sorted( Instances & instances, int W, OutIt out ) {
    typedef typename boost::range_value< Instances >::type w_mer_eval_t;
    boost::optional< w_mer_eval_t > last;
    BOOST_FOREACH( const w_mer_eval_t & instance, instances ) {

        // check sorted order
        BOOST_ASSERT( ! last || last->global_pos <= instance.global_pos );

        // is it the first one?
        if( ! last ) {
            last = instance;
        }
        // does it overlap with the last?
        else if( instance.global_pos < last->global_pos + typename w_mer_eval_t::global_pos_t( W ) ) {

            // keep the better one
            if( instance.Z > last->Z ) {
                last = instance;
            }

        } else {
            // otherwise output and update
            *out = *last;
            ++out;
            last = instance;
        }
    }
    if( last ) {
        *out = *last;
        ++out;
    }
}

/** Erase the overlapping positions from the vector. */
template< typename Spec >
void
erase_overlapping_from_score_sorted( typename w_mer_evaluation< Spec >::vec & instances, int W, int max_size ) {
    typename w_mer_evaluation< Spec >::vec::iterator i = instances.begin();
    while( instances.end() != i && i - instances.begin() != max_size ) {
        typename w_mer_evaluation< Spec >::vec::iterator j = instances.begin();
        for( ; j != i; ++j ) {
            if( wmer_overlap( i->global_pos, j->global_pos, W ) ) {
                instances.erase( i );
                break;
            }
        }
        if( i == j ) {
            ++i;
        }
    }
    instances.erase( i, instances.end() );
}



/** Output the non-overlapping positions from the range. */
template< typename InIt, typename OutIt >
inline
void
output_nonoverlapping_from_score_sorted( InIt begin, InIt end, int W, int max_size, OutIt out ) {
    InIt i = begin;
    int count = 0;
    while( end != i && count != max_size ) {
        InIt j = begin;
        for( ; j != i; ++j ) {
            if( wmer_overlap( i->global_pos, j->global_pos, W ) ) {
                break;
            }
        }
        if( i == j ) {
            *out = *i;
            ++out;
            ++count;
        }
        ++i;
    }
}



/**
 * Calculate the amount of overlap between the sequences.
 *
 * Sequences of W-mers must be non-overlapping and sorted by position.
 */
template< typename InIt1, typename InIt2 >
inline
unsigned
calculate_overlap( InIt1 b1, InIt1 e1, int W1, InIt2 b2, InIt2 e2, int W2 ) {
    using namespace boost::icl;
    typedef typename InIt1::value_type::global_pos_t global_pos_t;

    unsigned result = 0;
    // merge both the sequences
    while( b1 != e1 && b2 != e2 ) {
        const right_open_interval< global_pos_t > i1( b1->global_pos, b1->global_pos + W1 );
        const right_open_interval< global_pos_t > i2( b2->global_pos, b2->global_pos + W2 );

        result += length( i1 & i2 );

        // advance the the iterator that is most behind
        if( first( i1 ) < first( i2 ) ) {
            const global_pos_t previous_end = i1.upper();
            ++b1;
            if( b1 != e1 && previous_end > b1->global_pos ) {
                throw std::logic_error(
                    MS_MAKE_STRING(
                        "Sequence 1 not sorted by position or has overlaps: "
                        << previous_end << " > " << b1->global_pos
                    )
                );
            }
        } else {
            const global_pos_t previous_end = i2.upper();
            ++b2;
            if( b2 != e2 && previous_end > b2->global_pos ) {
                throw std::logic_error(
                    MS_MAKE_STRING(
                        "Sequence 2 not sorted by position or has overlaps: "
                        << previous_end << " > " << b2->global_pos
                    )
                );
            }
        }
    }
    return result;
}




/** Do instances overlap? */
GCC_DIAG_OFF(GCC_UNINIT_WARNING);
template< typename Instances >
bool
do_instances_overlap( Instances & instances, int W, bool already_sorted_by_position = false ) {
    typedef typename boost::range_value< Instances >::type w_mer_eval_t;
    if( ! already_sorted_by_position ) {
        sort_instances_by_position( instances );
    }
    boost::optional< w_mer_eval_t > last;
    BOOST_FOREACH( const w_mer_eval_t & instance, instances ) {
        BOOST_ASSERT( last->global_pos >= instance.global_pos );
        if( last && instance.global_pos < last->global_pos + typename w_mer_eval_t::global_pos_t( W ) ) {
            return true;
        }
        last = instance;
    }
    return false;
}
GCC_DIAG_ON(GCC_UNINIT_WARNING);




} // namespace steme

template< typename OS, typename Spec >
OS &
operator<<( OS & os, const steme::w_mer_evaluation< Spec > & e ) {
    return os << "Z=" << e.Z << "; pos=" << e.global_pos << ( e.rev_comp ? "-" : "+" );
}

#endif /* STEME_JR_6SEP2011_W_MER_EVALUATION_H_ */
