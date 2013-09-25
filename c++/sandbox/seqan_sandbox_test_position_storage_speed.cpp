/**
 * Copyright John Reid 2011
 *
 * \file
 * Code to test how efficient it is to store positions in the sequences as (seq id, offset) pairs or as global positions.
 */

#include <myrrh/seqan_boost_range_adaptors.h>

#include <seqan/index.h>
#include <seqan/sequence.h>

#include <boost/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>

#include <vector>
#include <stdexcept>

#include <seqan/version.h>
#if SEQAN_VERSION_MAJOR > 1 || SEQAN_VERSION_MINOR > 3
# define SEQAN_PACK Pack
# define SEQAN_BITPACKED BitPacked
#else
# define SEQAN_PACK Compressed
# define SEQAN_BITPACKED BitCompressed
#endif


using namespace seqan;

typedef String< Dna5 >                                         string_t;              /**< A string of Dna5. */
typedef StringSet< string_t >                                  string_set_t;          /**< StringSet type. */

/**
 * Read a FASTA file into a string set.
 */
void
read_fasta( const char * filename, string_set_t & sequences )
{
	::std::ifstream f( filename );
	if( ! f ) {
		std::cerr << "Could not open FASTA file: " << filename << "\n";
	} else {
		String< char > meta;
		string_t str;
		while( f ) {
			readMeta( f, meta, Fasta() );
			read( f, str, Fasta() );
			appendValue( sequences, str );
		}
	}
	f.close();
}


template< typename TIndex, typename Derived >
struct position_storage {

	typedef TIndex                                         index_t;
	typedef typename Host< index_t >::Type                 text_t;                ///< The type of the text in the index.
	typedef typename StringSetLimits< text_t >::Type       string_set_limits_t;   ///< Used to convert between local and global positions.
	typedef typename Iterator< TIndex, TopDown<> >::Type   top_down_it;           ///< A iterator over the index type.
	typedef typename Fibre< TIndex, EsaSA >::Type          sa_fibre_t;            ///< The type of the SA fibre.
	typedef typename Infix< sa_fibre_t const >::Type       occurrences_t;         ///< The type of a list of occurrences.
	typedef typename Size< occurrences_t >::Type           occ_size_t;            ///< The size type for a list of occurrences.
	typedef typename SAValue< TIndex >::Type               occurrence_t;          ///< The type of an occurrence.
	typedef size_t                                         seq_index_t;
	typedef size_t                                         offset_t;
	typedef typename Value< string_set_limits_t >::Type    global_pos_t;          ///< Type of a global position.
	typedef Pair< seq_index_t, offset_t >                  local_pos_t;

	/// Visit a node in the suffix tree (array).
	void
	visit( top_down_it it ) {

		if( repLength( it ) >= 8 ) {

			/// store positions if deeper than some depth
			//occurrences_t occs = seqan::getOccurrences( it );
			const occurrences_t & occs = getOccurrences( it );
			const occ_size_t num_occurrences = length( occs );
			for( occ_size_t i = 0; num_occurrences != i; ++i ) {
				static_cast< Derived * >( this )->store_position( occs[ i ] );
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


/// Stores positions as they are stored in the index
template< typename TIndex >
struct occurrence_storage : position_storage< TIndex, occurrence_storage< TIndex > > {
	typedef position_storage< TIndex, occurrence_storage< TIndex > > base;
	typedef std::vector< typename base::occurrence_t > occurrences_vec;

	occurrences_vec occurrences;

	occurrence_storage( TIndex & index ) { }

	inline
	void
	store_position( const typename base::occurrence_t & occ ) {
		occurrences.push_back( occ );
	}
};



/// Stores positions as (seq_id, offset) pairs
template< typename TIndex >
struct local_pos_storage : position_storage< TIndex, local_pos_storage< TIndex > > {
	typedef position_storage< TIndex, local_pos_storage< TIndex > > base;
	typedef std::vector< typename base::local_pos_t > occurrences_vec;

	occurrences_vec occurrences;

	local_pos_storage( TIndex & index ) { }

	inline
	void
	store_position( const typename base::occurrence_t & occ ) {
		occurrences.push_back( typename base::local_pos_t( getSeqNo( occ ), getSeqOffset( occ ) ) );
	}
};



/// Stores positions as global positions
template< typename TIndex >
struct global_pos_storage : position_storage< TIndex, global_pos_storage< TIndex > > {
	typedef position_storage< TIndex, global_pos_storage< TIndex > > base;
	typedef std::vector< typename base::global_pos_t > occurrences_vec;

	occurrences_vec occurrences;
	typename base::string_set_limits_t & string_set_limits;

	global_pos_storage( TIndex & index ) : string_set_limits( stringSetLimits( getFibre( index, EsaText() ) ) ) { }

	inline
	void
	store_position( const typename base::occurrence_t & occ ) {
		occurrences.push_back( posGlobalize( occ, string_set_limits ) );
	}
};



template< typename Storage >
boost::tuple< size_t, double >
time_storage( typename Storage::base::index_t & index ) {
	Storage storage( index );
	storage.occurrences.reserve( 13214001 );
	boost::timer timer;
	storage.visit( typename Storage::top_down_it( index ) );
	const double elapsed = timer.elapsed();
	return boost::make_tuple( storage.occurrences.size(), elapsed );
}


template< typename WMer >
int
dummy_wmer_fn( const WMer & wmer ) {
	int result = 0;
	typedef typename Iterator< WMer >::Type it_t;
	it_t e = end( wmer );
	for( it_t i = begin( wmer ); e != i; ++i ) {
		result += ordValue( *i );
	}
	return result;
}

struct local { };
struct compressed { };
struct global { };

typedef Index< string_set_t, IndexEsa< local > > local_index;
typedef Index< string_set_t, IndexEsa< compressed > > compressed_index;
typedef Index< string_set_t, IndexEsa< global > > global_index;

//struct local_index : Index< string_set_t > {
//	local_index( string_set_t & sequences ) : Index< string_set_t >( sequences ) { }
//};
//
//struct compressed_index : Index< string_set_t > {
//	compressed_index( string_set_t & sequences ) : Index< string_set_t >( sequences ) { }
//};
//
//struct global_index : Index< string_set_t > {
//	global_index( string_set_t & sequences ) : Index< string_set_t >( sequences ) { }
//};




namespace seqan {

// specialise for our indexes
template< typename TText >
struct SAValue< Index< TText, IndexEsa< compressed > > >
{
	typedef Pair< unsigned, unsigned, SEQAN_BITPACKED< 17, 13 > > Type;
};

// specialise for our indexes
template< typename TText >
struct SAValue< Index< TText, IndexEsa< global > > >
{
	typedef unsigned Type;
};

}


void
concatenate_sequences( string_set_t & sequences, string_t & concatenated_result ) {
	BOOST_FOREACH( string_t & seq, sequences ) {
		append( concatenated_result, seq );
	}
}


template< typename TIndex >
void
time_index_type( string_set_t & sequences, const char * name ) {

	std::cout << "Building " << name << " index.\n";
	std::cout << "sizeof( SAValue< " << name << " > ) = " << sizeof( typename SAValue< TIndex >::Type ) << "\n";
	boost::timer timer;
	TIndex index( sequences );
	// descend index once to make sure no delay due to initialisation
	time_storage< occurrence_storage< TIndex > >( index );
	std::cout << "Took " << timer.elapsed() << " seconds to build index.\n";

	std::cout << "Concatenating sequences\n";
	string_t concatenated;
	timer.restart();
	concatenate_sequences( sequences, concatenated );
	std::cout << "Took " << timer.elapsed() << " seconds to concatenate sequences.\n";

	//
	// Check global position storage
	//
	{
		size_t num_occurrences;
		double elapsed;
		boost::tie( num_occurrences, elapsed ) = time_storage< global_pos_storage< TIndex > >( index );
		std::cout
			<< "Global position storage with " << name << " index stored " << num_occurrences << " occurrences in "
			<< elapsed << " seconds.\n"
			;
	}

	//
	// Check local position storage
	//
	{
		size_t num_occurrences;
		double elapsed;
		boost::tie( num_occurrences, elapsed ) = time_storage< local_pos_storage< TIndex > >( index );
		std::cout
			<< "Local position storage with " << name << " index stored " << num_occurrences << " occurrences in "
			<< elapsed << " seconds.\n"
			;
	}

	//
	// Check occurrence storage
	//
	{
		size_t num_occurrences;
		double elapsed;
		boost::tie( num_occurrences, elapsed ) = time_storage< occurrence_storage< TIndex > >( index );
		std::cout
			<< "Occurrence storage with " << name << " index stored " << num_occurrences << " occurrences in "
			<< elapsed << " seconds.\n"
			;
	}


	//
	// Test access to W-mers.
	//
	const size_t W = 8;
	global_pos_storage< TIndex > storage( index );
	storage.occurrences.reserve( 13214001 );
	typename Iterator< TIndex, TopDown<> >::Type it2( index );
	typename global_pos_storage< TIndex >::top_down_it it( index );
	storage.visit( typename global_pos_storage< TIndex >::top_down_it( index ) );
	timer.restart();
	int total = 0;
	BOOST_FOREACH( typename global_pos_storage< TIndex >::global_pos_t pos, storage.occurrences ) {
		const typename global_pos_storage< TIndex >::seq_index_t seq = getSeqNo( pos, storage.string_set_limits );
		const typename global_pos_storage< TIndex >::offset_t offset = getSeqOffset( pos, storage.string_set_limits );
		total += dummy_wmer_fn( infix( getValueById( indexText( index ),  seq ), offset, offset + W ) );
	}
	std::cout << "Getting W-mers by local positions took " << timer.elapsed() << " seconds. Total was " << total << "\n";
	timer.restart();
	total = 0;
	BOOST_FOREACH( typename global_pos_storage< TIndex >::global_pos_t pos, storage.occurrences ) {
		total += dummy_wmer_fn( infix( indexRawText( index ), pos, pos+W ) );
	}
	std::cout << "Getting W-mers by global positions took " << timer.elapsed() << " seconds. Total was " << total << "\n";
	timer.restart();
	total = 0;
	BOOST_FOREACH( typename global_pos_storage< TIndex >::global_pos_t pos, storage.occurrences ) {
		total += dummy_wmer_fn( infix( concatenated, pos, pos+W ) );
	}
	std::cout << "Getting W-mers by global concatenated positions took " << timer.elapsed() << " seconds. Total was " << total << "\n";

}



int
main( int argc, char * argv[] ) {

	if( argc < 2 ) {
		std::cerr << "USAGE: " << argv[0] << " <fasta file>\n";
		return -1;
	} else {
		string_set_t sequences;
		boost::timer timer;
		read_fasta( argv[1], sequences );
		Size< string_t >::Type max_length = 0;
		Size< string_t >::Type total_bases = 0;
		for( Size< string_set_t >::Type i = 0; length( sequences ) != i; ++i ) {
			Size< string_t >::Type len = length( value( sequences, i ) );
			total_bases += len;
			max_length = std::max( max_length, len );
		}
		std::cout << "Read " << total_bases << " base pairs from " << length( sequences ) << " sequences from " << argv[ 1 ] << " in " << timer.elapsed() << " seconds\n";
		std::cout << "Longest sequence has " << max_length << " base pairs\n";

		time_index_type< local_index >( sequences, "local" );
		time_index_type< compressed_index >( sequences, "compressed" );
		//time_index_type< global_index >( sequences, "global" );

		return 0;
	}
}
