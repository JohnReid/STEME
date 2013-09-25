/**
 * Copyright John Reid 2011, 2012
 */

#include <steme/data.h>

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <boost/timer.hpp>

#ifdef USING_GOOGLE_PROFILER
#include <google/profiler.h>
#endif //USING_GOOGLE_PROFILER

#include <string>
#include <iostream>
#include <limits>

#include <seqan/version.h>
#if SEQAN_VERSION_MAJOR > 1 || SEQAN_VERSION_MINOR > 3
# define SEQAN_PACK Pack
# define SEQAN_BITPACKED BitPacked
#else
# define SEQAN_PACK Compressed
# define SEQAN_BITPACKED BitCompressed
#endif

using namespace seqan;
using namespace std;


typedef Dna5 dna_t;
typedef String< dna_t > string_t;
typedef StringSet< string_t > string_set_t;
typedef IndexEsa<> index_impl_t;
//typedef IndexWotd<> index_impl_t;
typedef Index< string_set_t, index_impl_t > index_t;
typedef Iterator< index_t, TopDown<> >::Type top_down_it;
typedef Iterator< index_t, TopDown< ParentLinks<> > >::Type top_down_history_it;


//
// Specialise STEME types
//
namespace steme {
/// Meta-function to choose type.
template<>
struct index_tag_meta< default_spec > {
    typedef index_impl_t type; ///< Selector for which seqan index implementation.
};
} // namespace steme

//
// Specialise various types to reduce memory requirements for ESA
//
namespace seqan
{
    template <>
    struct Fibre< index_t, FibreSA >
    {
        typedef unsigned int seq_index_t;
        typedef unsigned int offset_t;

        /**
         * Works for up to numeric_limits< seq_index_t >::max() sequences
         * of maximum length numeric_limits< offset_t >::max() base pairs
         */
        typedef String<
                        Pair< seq_index_t, offset_t, SEQAN_PACK >,
                        DefaultIndexStringSpec< index_t >::Type
                       > Type;

        // Use a mmapped string
        // typedef String< Pair<unsigned char, unsigned int, SEQAN_PACK>, MMap<> >   Type;
    };


    template <>
    struct Fibre< index_t, FibreLcp >
    {
        typedef String< unsigned int, DefaultIndexStringSpec< index_t>::Type >   Type;
    };


    template <>
    struct Fibre< index_t, FibreChildtab >
    {
        typedef String< unsigned int, DefaultIndexStringSpec< index_t>::Type >   Type;
    };
}



size_t
load_using_mmap( const char * filename, string_set_t & seqs ) {

    MultiSeqFile multi_seq_file;
    size_t num_bases = 0;

    // First we associate our sequence file with the memory mapped string underlying the ConcatDirect StringSet using open.
    open(multi_seq_file.concat, filename, OPEN_RDONLY);

    // Next we guess the file format of the single concatenation string and store the result in a AutoSeqFormat object,
    // which is used subsequently to select the right import function. split expects a ConcatDirect StringSet and divides
    // the underlying string into sequence fragments separated by a file format specific delimiter.
    AutoSeqFormat format;
    guessFormat(multi_seq_file.concat, format);
    split(multi_seq_file, format);

    // After calling split the multiSeqFile StringSet represents the sequence fragments and can be used to reserve memory
    // for the StringSets that store sequences and ids.
    unsigned seqCount = length(multi_seq_file);
    StringSet< CharString > seqIDs;
    reserve(seqs, seqCount, Exact());
    reserve(seqIDs, seqCount, Exact());

    // The main loop iterates over each sequence fragment and uses the functions assignSeq and assignSeqId to
    // extract sequence data and id.
    string_t seq;
    CharString id;
    for (unsigned i = 0; i < seqCount; ++i)
    {
        assignSeq(seq, multi_seq_file[i], format);    // read sequence
        assignSeqId(id, multi_seq_file[i], format);   // read sequence id

        // we use reserve and append, as assign is not supported
        // by StringSet<..., Owner<ConcatDirect<> > >
        appendValue(seqs, seq, Generous());
        appendValue(seqIDs, id, Generous());
        num_bases += length( seq );
    }
    return num_bases;
}




/**
 * Read a FASTA file into a string set.
 */
size_t
load_fasta( const char * filename, string_set_t & sequences )
{
    size_t num_bases = 0;
    ::ifstream f( filename );
    if( ! f ) {
        cerr << "Could not open FASTA file: " << filename << "\n";
    } else {
        String< char > meta;
        string_t str;
        while( f ) {
            readMeta( f, meta, Fasta() );
            read( f, str, Fasta() );
            appendValue( sequences, str );
            num_bases += length( str );
        }
    }
    f.close();
    return num_bases;
}




/**
 * Visit a node in the suffix tree (array).
 */
template< typename It >
void
visit( It it ) {
	if( goDown( it ) ) {
		visit( it );
		while( goRight( it ) ) {
			visit( it );
		}
	}
}


void wait_for_return() {
	char input;
	cout << "Press <return> to continue...\n";
	do {
		input = cin.get();
	} while( '\n' != input );
}

#pragma GCC diagnostic ignored "-Wunused-result"

void
show_ps_info() {
	system( "ps -C seqan_sandbox_build_index -o pid,sz,vsz,rss,args" );
}



int
main( int argc, char * argv[] ) {

	if( argc < 2 ) {

		cerr << "USAGE: " << argv[0] << " <fasta-file> [profile-filename]\n";
		return -1;

	} else {

#ifdef USING_GOOGLE_PROFILER
	    if( argc > 2 ) { ProfilerStart( argv[ 2 ] ); }
#endif //USING_GOOGLE_PROFILER

	    typedef Value< Fibre< index_t, FibreSA >::Type >::Type fibre_sa_t;
        typedef Value< fibre_sa_t, 1 >::Type seq_index_t;
        typedef Value< fibre_sa_t, 2 >::Type offset_t;
        const seq_index_t max_num_seqs   = numeric_limits< seq_index_t >::max();
        const offset_t    max_seq_length = numeric_limits< offset_t    >::max();
        cout << "Maximum # sequences = " << max_num_seqs << endl;
        cout << "Maximum sequence length = " << max_seq_length << endl;

	    show_ps_info();
	    boost::timer timer;

	    timer.restart();
	    string_set_t sequences;
	    const size_t num_bases = load_using_mmap( argv[ 1 ], sequences );
	    cout
	        << "\nRead " << num_bases
	        << " bases from " << length( sequences )
	        << " sequences in: " << argv[ 1 ]
	        << " in " << timer.elapsed()
	        << " seconds.\n";
	    show_ps_info();

	    // check index can hold sequences
	    typedef Size< string_set_t >::Type seqs_size_t;
        if( length( sequences ) > seqs_size_t( max_num_seqs ) ) {
            throw std::length_error( "Too many sequences." );
        }
        for( seqs_size_t i = 0; i != length( sequences ); ++i ) {
            if( length( value( sequences, i ) ) > max_seq_length ) {
                throw std::length_error( "Sequence too long." );
            }
        }

	    timer.restart();
	    index_t index( sequences );
	    cout << "\nBuilt index in " << timer.elapsed() << " seconds.\n";
	    show_ps_info();

        // first time could be slower as index is built on-the-fly
	    timer.restart();
	    visit( top_down_it( index ) );
	    cout << "\nVisited index in " << timer.elapsed() << " seconds.\n";
	    show_ps_info();

	    // second time should be much quicker as index has been built
        timer.restart();
        visit( top_down_it( index ) );
        cout << "\nVisited index in " << timer.elapsed() << " seconds.\n";
        show_ps_info();

        // just very slow with a TopDownHistory iterator
//        timer.restart();
//        visit( top_down_history_it( index ) );
//        cout << "\nVisited index in " << timer.elapsed() << " seconds with TopDownHistory iterator.\n";
//        show_ps_info();

	    timer.restart();
	    size_t max_W = 14;
	    steme::data< steme::default_spec > data( index, max_W );
	    cout << "\nBuilt data in " << timer.elapsed() << " seconds.\n";
	    show_ps_info();

#ifdef USING_GOOGLE_PROFILER
        if( argc > 2 ) { ProfilerStop(); }
#endif //USING_GOOGLE_PROFILER

		return 0;
	}
}
