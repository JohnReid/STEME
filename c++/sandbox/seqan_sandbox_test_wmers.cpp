/**
 * Copyright John Reid 2011
 *
 * \file
 * Code to test how to retrieve W-mers from their locations.
 */

#include <myrrh/seqan_boost_range_adaptors.h>

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io/guess_stream_format.h> // had to add these 2 lines for it to compile with latest
#include <seqan/seq_io/read_fasta_fastq.h>    // svn trunk: rev12932

#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <fstream>
#include <stdexcept>

using namespace seqan;
using namespace std;

typedef String< Dna5 >                        string_t;              ///< A string of Dna5.
typedef StringSet< string_t >                 string_set_t;          ///< StringSet type.
typedef StringSet<CharString>                 id_set_t;              ///< Holds IDs of strings.
typedef Index< string_set_t, IndexEsa<> >     index_t;               ///< Enhanced suffix array.
typedef Host< index_t >::Type                 text_t;                ///< The type of the text in the index.
typedef Fibre< index_t, EsaRawText >::Type    raw_text_t;            ///< The type of the raw text in the index.
typedef StringSetLimits< text_t >::Type       string_set_limits_t;   ///< Used to convert between local and global positions.
typedef Iterator< index_t, TopDown<> >::Type  top_down_it;           ///< A iterator over the index type.
typedef Fibre< index_t, EsaSA >::Type         sa_fibre_t;            ///< The type of the SA fibre.
typedef Infix< sa_fibre_t const >::Type       occurrences_t;         ///< The type of a list of occurrences.
typedef Size< occurrences_t >::Type           occ_size_t;            ///< The size type for a list of occurrences.
typedef SAValue< index_t >::Type              occurrence_t;          ///< The type of an occurrence.
typedef size_t                                seq_index_t;
typedef size_t                                offset_t;
typedef Value< string_set_limits_t >::Type    global_pos_t;          ///< Type of a global position.
typedef Pair< seq_index_t, offset_t >         local_pos_t;

/**
 * Read a FASTA file into a string set.
 */
void
read_fasta( const char * filename, id_set_t & ids, string_set_t & sequences )
{
	// Open file and create RecordReader.
	ifstream fasta( filename, std::ios_base::in | std::ios_base::binary );
	if(! fasta.good() )
		throw std::logic_error( "Could not open FASTA file." );
	RecordReader< ifstream, SinglePass<> > reader( fasta );
	// Define variables for storing the sequences and sequence ids.
	if( read2( ids, sequences, reader, Fasta() ) != 0 ) {
		throw std::logic_error( "ERROR reading FASTA." );
	}

}

void
print_seq_info( const char * filename, string_set_t & sequences ) {
	Size< string_t >::Type max_length = 0;
	Size< string_t >::Type total_bases = 0;
	for( Size< string_set_t >::Type i = 0; length( sequences ) != i; ++i ) {
		Size< string_t >::Type len = length( value( sequences, i ) );
		total_bases += len;
		max_length = std::max( max_length, len );
	}
	std::cout << "Read " << total_bases << " base pairs from " << length( sequences ) << " sequences from " << filename << "\n";
	std::cout << "Longest sequence has " << max_length << " base pairs\n";
}



struct visitor {

	vector< occurrence_t > occurrences;

	/// Visit a node in the suffix tree (array).
	void
	visit( top_down_it it ) {

		if( repLength( it ) >= 8 ) {

			/// store positions if deeper than some depth
			//occurrences_t occs = seqan::getOccurrences( it );
			const occurrences_t & occs = getOccurrences( it );
			const occ_size_t num_occurrences = length( occs );
			for( occ_size_t i = 0; num_occurrences != i; ++i ) {
				occurrences.push_back( occs[ i ] );
			}

		} else {

			/// Visit rest of tree below node
			if( goDown( it ) ) {
				visit( it );
				while( goRight( it ) ) {
					visit( it );
				}
			}
		}
	}
};



/**
 * Functor that complements seqan::Dna5.
 */
struct ComplementDna5 : std::unary_function< seqan::Dna5, seqan::Dna5 > {

	/** Returns the complement of its argument. */
	seqan::Dna5 operator()( seqan::Dna5 b ) const {
		seqan::Dna5 result;
		switch( seqan::ordValue( b ) ) {
		case 0: result.value = 3; break;
		case 1: result.value = 2; break;
		case 2: result.value = 1; break;
		case 3: result.value = 0; break;
		case 4: result.value = 4; break;
		default: throw std::invalid_argument( "Unknown Dna5 base." );
		}
		return result;
	}
};


/**
 * Transforms an iterator into a reverse complement iterator.
 */
template< typename It >
boost::transform_iterator< ComplementDna5, boost::reverse_iterator< It > >
make_rev_comp_iterator( It it ) {
	return
		boost::make_transform_iterator(
			boost::make_reverse_iterator( it ),
			ComplementDna5()
		);
}



/**
 * Transforms a boost.range into its reverse complement.
 */
template< typename Range >
struct rev_comp_range {
	typedef
		boost::iterator_range<
			boost::transform_iterator<
				ComplementDna5,
				boost::reverse_iterator<
					typename boost::range_iterator< const Range >::type
				>
			>
		>
		Type
		;

};


template< typename Range >
typename rev_comp_range< Range >::Type
make_rev_comp_range( const Range & range ) {
	return
		boost::make_iterator_range(
			make_rev_comp_iterator( boost::end  ( range ) ),
			make_rev_comp_iterator( boost::begin( range ) )
		);
}



#define MS_DEBUG_STRING_RANGE( var_name, range ) \
	std::string var_name; \
	std::copy( boost::begin( range ), boost::end( range ), std::back_insert_iterator< std::string >( var_name ) )




int
main( int argc, char * argv[] ) {

	if( argc < 2 ) {
		std::cerr << "USAGE: " << argv[0] << " <fasta file>\n";
		return -1;
	} else {

		//
		// Load sequences
		//
		string_set_t sequences;
		id_set_t ids;
		boost::timer timer;
		read_fasta( argv[1], ids, sequences );
		std::cout << "Took " << timer.elapsed() << " seconds to read FASTA file: " << argv[1] << "\n";
		print_seq_info( argv[1], sequences );


		//
		// Build index and store W-mer locations
		//
		visitor v;
		index_t index( sequences );
		v.visit( top_down_it( index ) );
		std::cout << "Iterated over " << v.occurrences.size() << " W-mers\n";


		//
		// Test methods to retrieve W-mers
		//
		string_set_limits_t & string_set_limits = stringSetLimits( getFibre( index, EsaText() ) );
		BOOST_FOREACH( const occurrence_t & occ, v.occurrences ) {
			const seq_index_t seq = getSeqNo( occ, string_set_limits );
			const offset_t offset = getSeqOffset( occ, string_set_limits );
			const int global_pos = posGlobalize( occ, string_set_limits );

			typedef Infix< string_t >::Type local_infix_t;
			const local_infix_t local_infix = infix( getValueById( indexText( index ),  seq ), offset, offset + 8 );
			const string_t local_str = local_infix;
			const rev_comp_range< local_infix_t >::Type local_rev_comp_infix = make_rev_comp_range( local_infix );
			BOOST_ASSERT( 8 == boost::size( local_infix ) );
			BOOST_ASSERT( 8 == boost::size( local_rev_comp_infix ) );
			MS_DEBUG_STRING_RANGE( local_rev_comp, make_rev_comp_range( local_infix ) );

			typedef Infix< raw_text_t >::Type global_infix_t;
			BOOST_ASSERT( global_pos + 8 <= int( length( indexRawText( index ) ) ) );
			const global_infix_t global_infix = infix( indexRawText( index ), global_pos, global_pos + 8 );
			const string_t global_str = global_infix;
			const rev_comp_range< Infix< raw_text_t >::Type >::Type global_rev_comp_infix = make_rev_comp_range( global_infix );
			boost::range_iterator< global_infix_t >::type begin = boost::begin( global_infix );
			boost::range_iterator< global_infix_t >::type end = boost::end( global_infix );
			BOOST_ASSERT( end > begin );
			BOOST_ASSERT( 8 == length( global_infix ) );
			BOOST_ASSERT( 8 == boost::size( global_infix ) );
			BOOST_ASSERT( 8 == boost::size( global_rev_comp_infix ) );
			MS_DEBUG_STRING_RANGE( global_rev_comp, global_rev_comp_infix );

			if( local_str != global_str ) {
				throw std::logic_error( "W-mers do not match." );
			}

			if( local_rev_comp != global_rev_comp ) {
				throw std::logic_error( "Reverse-complement W-mers do not match." );
			}
		}

		return 0;
	}
}


