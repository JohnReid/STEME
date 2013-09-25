/**
 * Copyright John Reid 2011
 *
 * \file
 * Code to test how to descend a suffix tree in a preferred order.
 */

#include <myrrh/seqan_boost_range_adaptors.h>

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io/guess_stream_format.h> // had to add these 2 lines for it to compile with latest
#include <seqan/seq_io/read_fasta_fastq.h>    // svn trunk: rev12932

#include <boost/timer.hpp>
#include <boost/foreach.hpp>

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



template< typename Derived >
struct visitor_base {
	vector< occurrence_t > occurrences;
	size_t W;

	visitor_base( size_t W ) : W( W ) { }

	void
	store_occurrences( top_down_it it ) {
		const occurrences_t & occs = getOccurrences( it );
		const occ_size_t num_occurrences = length( occs );
		for( occ_size_t i = 0; num_occurrences != i; ++i ) {
			occurrences.push_back( occs[ i ] );
		}
	}

	/// Visit a node in the suffix tree (array).
	void
	visit( top_down_it it ) {

		if( repLength( it ) >= W ) {

			/// store positions if deeper than some depth
			//occurrences_t occs = seqan::getOccurrences( it );
			store_occurrences( it );

		} else {

			static_cast< Derived * >( this )->descend( it );
		}
	}
};


/// Class that does a straightforward lexicographic descent of the tree.
struct visitor : visitor_base< visitor > {

	visitor( size_t W ) : visitor_base< visitor >( W ) { }

	void descend( top_down_it it ) {
		/// Visit rest of tree below node
		if( goDown( it ) ) {
			visit( it );
			while( goRight( it ) ) {
				visit( it );
			}
		}
	}
};


/// Class that does a descent of the tree in some preferred order at each depth.
struct visitor_preferred : visitor_base< visitor_preferred > {

	typedef std::vector< string_t > string_vec;

	string_vec preferred_bases;

	visitor_preferred( size_t W ) : visitor_base< visitor_preferred >( W ) { }

	/// Visit a node in the suffix tree (array).
	void
	descend( top_down_it it ) {

		/// Visit rest of tree below node in preferred order.
		BOOST_FOREACH( Dna5 base, preferred_bases[ repLength( it ) ] ) {
			//cout << "Descending: " << setw( W ) << representative( it ) << ": depth=" << rep_length << "; base=" << base << "\n";
			top_down_it i = it;
			if( goDown( i, base ) ) {
				//cout << "Went down from \"" << representative( it ) << "\" -> \"" << representative( i ) << "\"\n";
				visit( i );
			}
		}
	}
};


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
		// Build index and visit it to make sure it is built
		//
		index_t index( sequences );
		const size_t W = 8;
		visitor v( W );
		v.visit( top_down_it( index ) );


		//
		// Visit with preferred order method
		//
		const size_t num_to_show = 20;
		{
			boost::timer t;
			visitor_preferred v( W );
			v.preferred_bases.push_back( "CGTNA" );
			v.preferred_bases.push_back( "GTNAC" );
			v.preferred_bases.push_back( "TNACG" );
			v.preferred_bases.push_back( "ACGTN" );
			v.preferred_bases.push_back( "CGTNA" );
			v.preferred_bases.push_back( "GTNAC" );
			v.preferred_bases.push_back( "TNACG" );
			v.preferred_bases.push_back( "ACGTN" );
			v.visit( top_down_it( index ) );
			std::cout << "Iterated over " << v.occurrences.size() << " W-mers in " << t.elapsed() << " seconds\n";
			for( size_t i = 0; num_to_show != i; ++i ) {
				cout
					<< infixWithLength(
						getValue( indexText( index ), getSeqNo( v.occurrences[ i ] ) ),
						getSeqOffset( v.occurrences[ i ] ),
						W
					) << "\n"
					;
			}
		}



		//
		// Visit with normal method
		//
		{
			boost::timer t;
			visitor v( W );
			v.visit( top_down_it( index ) );
			std::cout << "Iterated over " << v.occurrences.size() << " W-mers in " << t.elapsed() << " seconds\n";
		}



	}
}
