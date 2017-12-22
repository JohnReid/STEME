/** Copyright John Reid 2011, 2012, 2013
 *
 * \file
 * \brief Tests branch and bound part of STEME algorithm.
 *
 */

#include <boost/assign/list_of.hpp>
#include <boost/test/utils/wrap_stringstream.hpp>
#if (BOOST_VERSION / 100) >= 159
# include <boost/test/included/unit_test.hpp>
#endif
#include <boost/test/floating_point_comparison.hpp>
#include <boost/timer.hpp>
#include <boost/filesystem.hpp>

#include <steme/data.h>
#include <steme/model.h>
#include <steme/descender.h>
#include <steme/seqan_types.h>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io/guess_stream_format.h> // had to add these 2 lines for it to compile with latest
#include <seqan/seq_io/read_fasta_fastq.h>    // svn trunk: rev12932

#include <fstream>

#define MAKE_STRING( x ) ( boost::wrap_stringstream().ref() << x ).str()

using namespace seqan;
using namespace std;
using namespace steme;
using boost::assign::list_of;
using namespace boost::filesystem;

typedef steme_seqan_types<> seqan_types_t;

typedef seqan_types_t::id_set_t                     id_set_t;
typedef seqan_types_t::string_set_t                 string_set_t;
typedef seqan_types_t::string_t                     string_t;
typedef seqan_types_t::text_t                       text_t;                ///< The type of the text in the index.
typedef seqan_types_t::raw_text_t                   raw_text_t;            ///< The type of the raw text in the index.
typedef seqan_types_t::string_set_limits_t          string_set_limits_t;   ///< Used to convert between local and global positions.
typedef seqan_types_t::occurrences_t                occurrences_t;         ///< The type of a list of occurrences.
typedef seqan_types_t::local_pos_t                  local_pos_t;

typedef model<>::bg_model_t                         bg_model_t;            ///< The background model type.
typedef model<>::bs_model_t                         bs_model_t;            ///< The binding site model type.


/**
 * Read a FASTA file into a string set.
 */
void
read_fasta(
    const path & filename,
    id_set_t & ids,
    string_set_t & sequences
) {
    // Open file and create RecordReader.
    std::ifstream fasta( filename.c_str(), std::ios_base::in | std::ios_base::binary );
    if(! fasta.good() )
        throw std::logic_error( "Could not open FASTA file." );
    RecordReader< std::ifstream, SinglePass<> > reader( fasta );
    // Define variables for storing the sequences and sequence ids.
    if( read2( ids, sequences, reader, Fasta() ) != 0 ) {
        throw std::logic_error( "ERROR reading FASTA." );
    }

}



void
print_seq_info( const path & filename, string_set_t & sequences ) {
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


unsigned
get_num_above_threshold( const std::vector< double > & calculated, double threshold ) {
    unsigned result = 0;
    BOOST_FOREACH( double x, calculated ) {
        result += ( x >= threshold );
    }
    return result;
}


void
throw_if_not_close( double x, double y, double tolerance, const char * msg ) {
#if (BOOST_VERSION / 100) >= 159
    // quick hack to get working with modern boost versions, will not throw error like old implementation
    BOOST_REQUIRE_CLOSE( x, y, tolerance );
#else
    using ::boost::test_tools::check_is_close;
    using ::boost::test_tools::fraction_tolerance;
    using ::boost::test_tools::percent_tolerance;
    if( ! check_is_close( x, y, fraction_tolerance( tolerance ) ) ) {
        throw std::logic_error(
            MYRRH_MAKE_STRING(
                msg << ": Not close, " << x << "!=" << y
                << ", difference=" << x-y << " (fraction=" << (x-y)/y << ")"
            ).c_str()
        );
    }
#endif
}



struct visitor : tree_descender< visitor > {

    double log_threshold;
    std::vector< double > positive_Z;
    std::vector< double > negative_Z;
    std::vector< double > positive_log_eval;
    std::vector< double > negative_log_eval;

    /** Constructor. */
    visitor(
        data<> &         _data,
        model<> &        _model,
        double          log_threshold
    )
    : tree_descender< visitor >( _data, _model )
    , log_threshold( log_threshold )
    , positive_Z( _data.N, -1. )
    , negative_Z( _data.N, -1. )
    , positive_log_eval( _data.N, 1. )
    , negative_log_eval( _data.N, 1. )
    { }

    unsigned
    num_above_threshold( double threshold ) const {
        return
            get_num_above_threshold( positive_Z, threshold )
            + get_num_above_threshold( negative_Z, threshold )
            ;
    }

    /**
     * \return The threshold.
     */
    double
    get_llr_threshold() {
        return log_threshold;
    }


    GCC_DIAG_OFF(uninitialized);
    template< bool RevComp >
    inline
    void
    handle_W_mer_occurrence_one_strand(
        int global_pos,
        double Z_strand,
        double g_n,
        double lambda,
        double log_lambda,
        double log_p_X_n_theta_BS,
        double log_p_background
    ) {
        const double model_calculated_Z = _model.calculate_Z( global_pos, RevComp );
        //throw_if_not_close( 1.51, 1.50, 0.005, "Test" );
        throw_if_not_close( model_calculated_Z, Z_strand, 1e-12, "Model calculated Z != tree descender calculated Z" );
        if( RevComp ) {
            negative_Z[ global_pos ] = Z_strand;
            negative_log_eval[ global_pos ] = log_p_X_n_theta_BS;
        } else {
            positive_Z[ global_pos ] = Z_strand;
            positive_log_eval[ global_pos ] = log_p_X_n_theta_BS;
        }
    }
    GCC_DIAG_ON(uninitialized);


    /// Handle a calculated Z for both strands of an occurrence of a W-mer
    inline
    void
    handle_W_mer_all_occurrences(
        top_down_it it,
        size_t num_occurrences,
        double Z_positive_total,
        double Z_negative_total
    ) {
    }
};

namespace steme {
const double log_quarter = std::log( 0.25 ); /**< log(.25) */
}

void
check_Z(
    data<> & _data,
    size_t W,
    size_t pos,
    double Z_threshold,
    double Z,
    double calculated_Z,
    bool rev_comp
) {
    if( Z >= Z_threshold && calculated_Z != Z ) {
        std::string w_mer;
        assign( w_mer, _data.get_W_mer( W, pos ) );
        throw std::logic_error(
            MAKE_STRING(
                "Missing a Z on " << ( rev_comp ? "negative" : "positive" ) << " strand at global position " << pos
                << "; W-mer=" << w_mer
                << ": Z=" << Z
                << "; calculated Z=" << calculated_Z
                << "; threshold=" << Z_threshold
            )
        );
    }
}


void
test_seed(
    data<> & _data,
    zero_order_frequencies & freqs_with_pseudo_counts,
    likelihoods_vec_vec_t & lls,
    string_t seed
) {
    std::cout << "Testing seed: " << seed << "\n";

    const size_t W = length( seed );
    bg_model_t bg( W, freqs_with_pseudo_counts, &lls );
    base_to_wmer_likelihoods_calculator( _data, *bg.base_LLs, bg.wmer_LLs );

    //
    // Build model
    //
    bs_model_t bs( Pssm( PssmStorage( W, 4 ) ) );
    bs.seed( seed );
    model<> _model( _data, bs, bg );


    //
    // Descend tree to find all Zs using a threshold of 0.
    //
    visitor all_Zs( _data, _model, -1e8 );
    all_Zs.descend_tree();
    BOOST_ASSERT( all_Zs.positive_Z.size() == all_Zs.negative_Z.size() );
    std::cout
        << "Have calculated " << all_Zs.num_above_threshold( -.5 )
        << " Z in total.\n";


    //
    // Descend tree for a range of thresholds and check all the Zs we were promised are there
    //
    const std::vector< double > thresholds = list_of
        ( -1.0e2 )
        ( -0.7e2 )
        ( -0.4e2 )
        ( -0.2e2 )
        ( -1.6e1 )
        ( -1.3e1 )
        ( -1.0e1 )
        ( -0.9e1 )
        ( -0.8e1 )
        ( -0.7e1 )
        ( -0.4e1 )
        ( -0.2e1 )
        ( -1.0e0 )
        ( -0.7e0 )
        ( -0.4e0 )
        ( -0.2e0 )
        (  0.0e0 )
        (  0.2e0 )
        (  0.4e0 )
        (  0.7e0 )
        (  1.0e0 )
        (  0.2e1 )
        (  0.4e1 )
        (  0.7e1 )
        (  1.0e1 )
        (  0.2e2 )
        ;
    BOOST_FOREACH( double threshold, thresholds ) {
        const double Z_threshold = lr_to_prob( std::exp( threshold ) );
        visitor v( _data, _model, threshold );
        v.descend_tree();
        const unsigned num_calculated = v.num_above_threshold( -.5 );
        const unsigned above_threshold = v.num_above_threshold( Z_threshold );
        std::cout
            << "Threshold=" << std::setw( 12 ) << Z_threshold
            << " : Have calculated Z for " << std::setw( 6 ) << num_calculated << " W-mers. "
            << std::setw( 6 ) << above_threshold << " were above threshold.\n";
        for( size_t pos = 0; all_Zs.positive_Z.size() != pos; ++pos ) {
            check_Z( _data, W, pos, Z_threshold, all_Zs.positive_Z[ pos ], v.positive_Z[ pos ], false );
            check_Z( _data, W, pos, Z_threshold, all_Zs.negative_Z[ pos ], v.negative_Z[ pos ], true  );
        }
    }
}


int
main( int argc, char * argv[] ) {


    //const char * fasta_filename = "/home/john/Data/MEIS1/meis1-peaks.fasta";
    path fasta_path;
    if( argc > 1 ) {
        fasta_path = path( argv[ 1 ] );
    } else {
        // find FASTA directory by going up directory hierarchy to STEME directory
        // and then down to FASTA directory
        path p = current_path();
        while( ! p.empty() && path( "STEME" ) != p.filename() ) {
            p = p.parent_path();
        }
        fasta_path = p / "python" / "stempy" / "test" / "fasta" / "T00759-tiny.fa";
    }

    std::vector< string_t > seeds = list_of
        ( "ACGT" )
        ( "TATCT" )
        ( "GTAGA" )
        ( "TCAATG" )
        ( "AATTCCGG" )
        ( "CCAGCTTCT" )
        ( "CGGTGACTGGTTTTCTGCTCTCTCTC" )
        ;

    //
    // Load sequences
    //
    string_set_t sequences;
    id_set_t ids;
    boost::timer timer;
    read_fasta( fasta_path, ids, sequences );
    std::cout << "Took " << timer.elapsed() << " seconds to read FASTA file: " << fasta_path << "\n";
    print_seq_info( fasta_path, sequences );


    //
    // Build index
    //
    seqan_types_t::index_t index( sequences );


    //
    // Build stem object
    //
    data<> _data( index );


    //
    // We want to test branch-and-bound when using position-specific priors,
    // so set up some artifical priors by planting some TFBSs of varying
    // strengths.
    //
    const Size< string_t >::Type psp_W = 14;
    const Size< string_t >::Type spacing = 37;
    for( Size< string_set_t >::Type seq = 0; length( sequences ) != seq; ++seq ) {
        const Size< string_t >::Type len = length( value( sequences, seq ) );
        for( Size< string_t >::Type pos = psp_W / 2; pos < len - psp_W; pos += spacing ) {
            const double Z = 1. - ( pos % 9 ) * .1; // some randomish strength
            _data._base_prior->adjust_for_binding_site( seq, pos, psp_W, Z );
        }
    }


    //
    // Build background
    //
    typedef complete_markov_model< 2, 4, double > markov_model_t;
    markov_model_t mm;
    zero_order_frequencies freqs = build_model_from_index( mm, index, 1. );
    zero_order_frequencies freqs_with_pseudo_counts = freqs.add_pseudo_counts( 1. );
    likelihoods_vec_vec_t lls;
    calculate_likelihoods( mm, _data.get_text(), lls );

    BOOST_FOREACH( string_t & seed, seeds ) {
        test_seed( _data, freqs_with_pseudo_counts, lls, seed );
    }

}

