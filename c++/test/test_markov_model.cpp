/**
 * Copyright John Reid 2011, 2012
 */

#define BOOST_TEST_MODULE markov_model test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/timer.hpp>

#include <steme/markov.h>
#include <steme/background_model.h>
#include <steme/seqan_types.h>
#include <vector>
#include <boost/assign/list_of.hpp>

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io/guess_stream_format.h> // had to add these 2 lines for it to compile with latest
#include <seqan/seq_io/read_fasta_fastq.h>    // svn trunk: rev12932


//#include <seqan/version.h>
//#if SEQAN_VERSION_MAJOR > 1 || SEQAN_VERSION_MINOR > 3
//# define SEQAN_READ read
//#else
//# define SEQAN_READ read2
//#endif


using namespace steme;
using namespace boost::assign;
using namespace std;
using namespace seqan;

typedef steme_seqan_types<> seqan_types_t;

BOOST_AUTO_TEST_CASE( multi_array_access_test )
{
	// our array
	typedef vector< size_t > vector_t;
	boost::multi_array< double, 4 > a( boost::extents[2][2][2][2] );
	vector< size_t > index = list_of(0)(1)(1)(0);
	a[0][1][1][0] = 7.;
	using steme::detail::access_multi_array_element;
	BOOST_CHECK_EQUAL( a[0][1][1][0], access_multi_array_element( boost::type< double >(), a, index.begin() ) );
	BOOST_CHECK_EQUAL( a[0][1][1][0], a( index ) );
}




BOOST_AUTO_TEST_CASE( markov_model_init_test )
{
	typedef complete_markov_model< 3, 4, double > model_t;

	{
		model_t model;
		BOOST_CHECK_EQUAL( model.mm.storage[3][1][2][0], 0. );

		// initialise counts with pseudo-count of 6.
		add_pseudo_counts( model, double( 6 ) );
		BOOST_CHECK_EQUAL( model.mm.storage[3][1][2][0], 6. );
		BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][0], 6. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][0], 6. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[0], 6. );

		// add pseudo-counts of 3.
		add_pseudo_counts( model, double( 3 ) );
		BOOST_CHECK_EQUAL( model.mm.storage[3][1][2][0], 9. );
	}

	{
		model_t model;
		// add counts from a sequence
		vector< size_t > seq = list_of(0)(1)(2)(3)(0);
		add_counts_to_complete( model, seq.begin(), seq.end() );
		BOOST_CHECK_EQUAL( model.mm.storage[0][1][2][3], 1. );
		BOOST_CHECK_EQUAL( model.mm.storage[3][1][2][3], 0. );
		BOOST_CHECK_EQUAL( model.mm.storage[0][0][2][3], 0. );
		BOOST_CHECK_EQUAL( model.mm.storage[0][1][1][3], 0. );
		BOOST_CHECK_EQUAL( model.mm.storage[0][1][2][2], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][3], 1. );
		BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[3][2][3], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][3][3], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][0], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][3], 1. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[1][3], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][0], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[0], 2. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[1], 1. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[2], 1. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[3], 1. );

		// normalise counts
		normalise_counts( model );
		BOOST_CHECK_EQUAL( model.mm.storage[0][1][2][3], 1. );
		BOOST_CHECK_MESSAGE( std::isnan( model.mm.storage[3][1][2][3] ), "Should be NaN" );
		BOOST_CHECK_MESSAGE( std::isnan( model.mm.storage[0][0][2][3] ), "Should be NaN" );
		BOOST_CHECK_MESSAGE( std::isnan( model.mm.storage[0][1][1][3] ), "Should be NaN" );
		BOOST_CHECK_EQUAL( model.mm.storage[0][1][2][2], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][3], 1. );
		BOOST_CHECK_MESSAGE( std::isnan( model.lower_orders.mm.storage[3][2][3] ), "Should be NaN" );
		BOOST_CHECK_MESSAGE( std::isnan( model.lower_orders.mm.storage[1][3][3] ), "Should be NaN" );
		BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][0], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][3], 1. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[1][3], 0. );
		BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][0], 0. );
		BOOST_CHECK_CLOSE( model.lower_orders.lower_orders.lower_orders.mm.storage[0], 2./5., 0.0001 );

		// convert to scale
		convert_to_scale( model );
		BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][3], 0. );
		BOOST_CHECK_CLOSE( model.lower_orders.lower_orders.lower_orders.mm.storage[0], std::log( 2./5. ), 0.0001 );
	}
}



BOOST_AUTO_TEST_CASE( markov_model_evaluation_test )
{
	// create and initialise model.
	typedef complete_markov_model< 3, 4, double > model_t;
	model_t model;
	vector< size_t > seq = list_of(0)(1)(2)(3)(0);
	initialise_from_sequence( model, 1., seq.begin(), seq.end() );

	{
		vector< size_t > eval_seq = list_of(0);
		vector< double > probs;
		model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
		BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
		BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
	}

	{
		vector< size_t > eval_seq = list_of(0)(1);
		vector< double > probs;
		model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
		BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
		BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
	}

	{
		vector< size_t > eval_seq = list_of(0)(1)(2);
		vector< double > probs;
		model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
		BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
		BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
	}

	{
		vector< size_t > eval_seq = list_of(0)(1)(2)(3);
		vector< double > probs;
		model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
		BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
		BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
	}

	{
		vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0);
		vector< double > probs;
		model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
		BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
		BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 2./5. ), 0.0001 );
	}

	{
		vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1);
		vector< double > probs;
		model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
		BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
		BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 2./5. ), 0.0001 );
		BOOST_CHECK_CLOSE( probs[5] - probs[4], std::log( 1./4. ), 0.0001 );
	}

    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
        BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[5] - probs[4], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[6] - probs[5], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[7] - probs[6], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is not large enough
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 1,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 1, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[5] - probs[4], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[6] - probs[5], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is not large enough
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 2,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 2, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[5] - probs[4], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is just right size
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 3,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 3, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is large enough
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 4,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 4, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is large enough
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 6,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 6, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
    }
}


/**
 * Read a FASTA file into a string set.
 */
void
read_fasta(
	const char * filename,
	seqan_types_t::id_set_t & ids,
	seqan_types_t::string_set_t & sequences
) {
	// Open file and create RecordReader.
	ifstream fasta( filename, std::ios_base::in | std::ios_base::binary );
	if(! fasta.good() )
		throw std::logic_error( "Could not open FASTA file." );
	RecordReader< ifstream, DoublePass<> > reader( fasta );
	// Define variables for storing the sequences and sequence ids.
	if( read2( ids, sequences, reader, Fasta() ) != 0 ) {
		throw std::logic_error( "ERROR reading FASTA." );
	}

}

void
print_seq_info(
	const char * filename,
	seqan_types_t::string_set_t & sequences
) {
	Size< seqan_types_t::string_t >::Type max_length = 0;
	Size< seqan_types_t::string_t >::Type total_bases = 0;
	for( Size< seqan_types_t::string_set_t >::Type i = 0; length( sequences ) != i; ++i ) {
		Size< seqan_types_t::string_t >::Type len = length( value( sequences, i ) );
		total_bases += len;
		max_length = std::max( max_length, len );
	}
	std::cout << "Read " << total_bases << " base pairs from " << length( sequences ) << " sequences from " << filename << "\n";
	std::cout << "Longest sequence has " << max_length << " base pairs\n";
}



BOOST_AUTO_TEST_CASE( markov_model_build_test )
{
	const char * fasta_filename = "../python/stempy/test/fasta/T00759-tiny.fa";

	//
	// Load sequences
	//
	seqan_types_t::string_set_t sequences;
	seqan_types_t::id_set_t ids;
	boost::timer timer;
	read_fasta( fasta_filename, ids, sequences );
	std::cout << "Took " << timer.elapsed() << " seconds to read FASTA file: " << fasta_filename << "\n";
	print_seq_info( fasta_filename, sequences );

	//
	// Build index
	//
	seqan_types_t::index_t index( sequences );

	//
	// create and initialise models
	//
	typedef complete_markov_model< 3, 4, double > model_t;
	model_t model_from_seqs;
	add_counts_from_seqs_to_complete( model_from_seqs, sequences );
	model_t model_from_index;
	add_counts_from_index_to_complete( model_from_index, index );
	BOOST_CHECK( model_from_seqs.mm.storage == model_from_index.mm.storage );
}

