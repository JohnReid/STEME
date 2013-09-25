/**
 * Copyright John Reid 2011, 2012
 *
 * \file Code to test speed of different find best W-mer implementations.
 */

// turn timer info on.
#define STEME_LOG_TIMING_INFO

#include <steme/find_best_w_mers.h>

#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/lambda/lambda.hpp>

#include <seqan/stream.h>
#include <seqan/index/index_sa_lss.h>
#include <seqan/seq_io/guess_stream_format.h> // had to add these 2 lines for it to compile with latest
#include <seqan/seq_io/read_fasta_fastq.h>    // svn trunk: rev12932

using namespace boost;
using namespace std;
using namespace steme;
using namespace seqan;


typedef string_meta<>::type string_t;
typedef string_set_meta<>:: type string_set_t;
typedef id_set_meta<>:: type id_set_t;
typedef index_meta<>::type index_t;
typedef Iterator< index_t, TopDown<> >::Type top_down_it;
typedef Fibre< index_t, EsaSA >::Type sa_fibre_t;
typedef Infix< sa_fibre_t const >::Type occurrences_t;
typedef Size< occurrences_t >::Type occ_size_t;
typedef SAValue< index_t >::Type occurrence_t;


namespace steme {
const double log_quarter = std::log( 0.25 ); /**< log(.25) */
}


/**
 * Read a FASTA file into a string set.
 */
void
read_fasta( const char * filename, id_set_t & ids, string_set_t & sequences )
{
    // Open file and create RecordReader.
    ifstream fasta( filename, ios_base::in | ios_base::binary );
    if(! fasta.good() )
        throw logic_error( "Could not open FASTA file." );
    RecordReader< ifstream, SinglePass<> > reader( fasta );
    // Define variables for storing the sequences and sequence ids.
    if( read2( ids, sequences, reader, Fasta() ) != 0 ) {
        throw logic_error( "ERROR reading FASTA." );
    }

}

void
print_seq_info( const char * filename, string_set_t & sequences ) {
    Size< string_t >::Type max_length = 0;
    Size< string_t >::Type total_bases = 0;
    for( Size< string_set_t >::Type i = 0; length( sequences ) != i; ++i ) {
        Size< string_t >::Type len = length( value( sequences, i ) );
        total_bases += len;
        max_length = max( max_length, len );
    }
    cout << "Read " << total_bases << " base pairs from " << length( sequences ) << " sequences from " << filename << "\n";
    cout << "Longest sequence has " << max_length << " base pairs\n";
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

template< typename It >
void
show_n( It begin, It end, unsigned n ) {
    for( unsigned i = 0; n != i && end != begin; ++i ) {
        cout << *begin++ << "\n";
    }
}


template< typename FindBestWMers >
void
time_find_best_w_mers( data<> & _data, model<> & _model, const string & name, size_t num_to_find ) {
    using namespace boost::adaptors;

    cout << "***************** " << name << " method ******************\n";
    boost::timer t;
    FindBestWMers finder( _data, _model, num_to_find );
    finder.descend_tree();
    cout << "Took " << t.elapsed() << " seconds to find best W-mers.\n";
    cout << "Found " << finder.best_w_mers->size() << " W-mers.\n";
    if( finder.has_overlapping() ) {
        cout << "!!!!!!!!!!!!!!!!!!!!!!! HAS OVERLAPS !!!!!!!!!!!!!!!!!!!!!\n";
    }
    unsigned num_to_show = 5;
    cout << "The " << num_to_show << " best W-mers:\n";
    show_n( finder.best_w_mers->begin(), finder.best_w_mers->end(), num_to_show );
    cout << "The " << num_to_show << " worst W-mers:\n";
    show_n( finder.best_w_mers->rbegin(), finder.best_w_mers->rend(), num_to_show );
}


int
main( int argc, char * argv[] ) {

    if( argc < 4 ) {

        cerr << "USAGE: " << argv[0] << " <fasta file> <seed> <number to find>\n";
        return -1;

    } else {

        const string fasta_file = argv[ 1 ];
        const string_t seed = argv[ 2 ];
        const size_t num_to_find = lexical_cast< size_t >( argv[ 3 ] );

        //
        // Load sequences
        //
        string_set_t sequences;
        id_set_t ids;
        boost::timer timer;
        read_fasta( fasta_file.c_str(), ids, sequences );
        cout << "Took " << timer.elapsed() << " seconds to read FASTA file: " << fasta_file << "\n";
        print_seq_info( fasta_file.c_str(), sequences );

        //
        // Build index
        //
        timer.restart();
        index_t index( sequences );
        cout << "Took " << timer.elapsed() << " seconds to build index.\n";

        //
        // Visit index to make sure it is built
        //
        timer.restart();
        const size_t W = length( seed );
        visitor v( W );
        v.visit( top_down_it( index ) );
        cout << "Took " << timer.elapsed() << " seconds to visit index.\n";


        //
        // Build a data object.
        //
        timer.restart();
        data<> _data( index, false, W );
        cout << "Took " << timer.elapsed() << " seconds to initialise data.\n";


        //
        // Build a Markov model of the sequences.
        //
        timer.restart();
        complete_markov_model< 2, 4, double > mm;
        zero_order_frequencies freqs = build_model_from_index( mm, index, 1. );
        cout << "Took " << timer.elapsed() << " seconds to build Markov model.\n";


        //
        // Build a background model.
        //
        timer.restart();
        likelihoods_vec_vec_t log_likelihoods;
        calculate_likelihoods( mm, _data.get_text(), log_likelihoods );
        background_model_meta<>::type bg( W, freqs, &log_likelihoods );
        base_to_wmer_likelihoods_calculator( _data, *bg.base_LLs, bg.wmer_LLs );

        cout << "Took " << timer.elapsed() << " seconds to build background model.\n";


        //
        // Build a binding model.
        //
        PssmStorage pssm_storage( W );
        Pssm pssm( pssm_storage );
        binding_model_meta<>::type bs( pssm );
        bs.seed( seed );


        //
        // Build a model.
        //
        model<> _model( _data, bs, bg );


        //
        // Find W-mers with sorted vector.
        //
        time_find_best_w_mers< find_best_w_mers< store_in_heap > >( _data, _model, "heap", num_to_find );
        time_find_best_w_mers< find_best_w_mers< store_in_sorted_vector > >( _data, _model, "sorted vector", num_to_find );
        time_find_best_w_mers< find_best_w_mers< store_all > >( _data, _model, "store all", num_to_find );
        //time_find_best_w_mers< find_best_w_mers< store_in_sorted_list > >( _data, _model, "sorted list", num_to_find );
        //time_find_best_w_mers< find_best_w_mers< store_using_multi_index > >( _data, _model, "multi-index", num_to_find );
        //time_find_best_w_mers< find_best_w_mers< store_in_set > >( _data, _model, "set", num_to_find );


        return 0;
    }
}

