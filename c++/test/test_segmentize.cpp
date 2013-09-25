/**
 * Copyright John Reid 2011
 */

#define BOOST_TEST_MODULE segmentize test
#include <boost/test/unit_test.hpp>
//#include <boost/test/floating_point_comparison.hpp>

#include <steme/background_model.h>
#include <vector>
#include <boost/assign/list_of.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/lambda/lambda.hpp>
#include <seqan/sequence.h>

using namespace steme;
using namespace boost;
using namespace boost::assign;
using namespace std;
using namespace seqan;


typedef vector< int > seq_t;
typedef seq_t::const_iterator iterator_t;
typedef boost::tuple< iterator_t, iterator_t > segment_t;
typedef vector< segment_t > segment_vec;

template< typename SegmentVec >
struct
segment_appender {
	SegmentVec & segments;

	segment_appender( SegmentVec & segments ) : segments( segments ) { }

	template< typename Iterator >
	void
	operator()( Iterator begin, Iterator end ) {
		segments.push_back( typename SegmentVec::value_type( begin, end ) );
	}
};

template< typename SegmentVec >
segment_appender< SegmentVec >
make_segment_appender( SegmentVec & segments ) {
	return segment_appender< SegmentVec >( segments );
}


void segment_seq( seq_t & seq, segment_vec & true_segments, segment_vec & false_segments ) {
	segmentize_sequence(
		seq.begin(),
		seq.end(),
		*boost::lambda::_1,
		make_segment_appender( true_segments ),
		make_segment_appender( false_segments )
	);

//	std::cout << "True segments\n";
//	BOOST_FOREACH( segment_t segment, true_segments ) {
//		std::cout << segment.get<0>()-seq.begin() << ":" << segment.get<1>()-seq.begin() << "\n";
//	}
//	std::cout << "False segments\n";
//	BOOST_FOREACH( segment_t segment, false_segments ) {
//		std::cout << segment.get<0>()-seq.begin() << ":" << segment.get<1>()-seq.begin() << "\n";
//	}
}


BOOST_AUTO_TEST_CASE( segmentize_1 )
{
	seq_t seq = list_of(0)(1)(1)(0);
	segment_vec true_segments;
	segment_vec false_segments;
	segment_seq( seq, true_segments, false_segments );

	BOOST_CHECK_EQUAL( int( true_segments.size() ), 1 );
	BOOST_CHECK( true_segments[0] == segment_t( seq.begin() + 1, seq.begin() + 3 ) );
	BOOST_CHECK_EQUAL( int( false_segments.size() ), 2 );
	BOOST_CHECK( false_segments[0] == segment_t( seq.begin() + 0, seq.begin() + 1 ) );
	BOOST_CHECK( false_segments[1] == segment_t( seq.begin() + 3, seq.begin() + 4 ) );
}


BOOST_AUTO_TEST_CASE( segmentize_2 )
{
	seq_t seq = list_of(1)(1)(1)(1);
	segment_vec true_segments;
	segment_vec false_segments;
	segment_seq( seq, true_segments, false_segments );

	BOOST_CHECK_EQUAL( int( true_segments.size() ), 1 );
	BOOST_CHECK( true_segments[0] == segment_t( seq.begin() + 0, seq.begin() + 4 ) );
	BOOST_CHECK_EQUAL( int( false_segments.size() ), 0 );
}


BOOST_AUTO_TEST_CASE( segmentize_3 )
{
	seq_t seq = list_of(0)(0)(0)(0);
	segment_vec true_segments;
	segment_vec false_segments;
	segment_seq( seq, true_segments, false_segments );

	BOOST_CHECK_EQUAL( int( false_segments.size() ), 1 );
	BOOST_CHECK( false_segments[0] == segment_t( seq.begin() + 0, seq.begin() + 4 ) );
	BOOST_CHECK_EQUAL( int( true_segments.size() ), 0 );
}

BOOST_AUTO_TEST_CASE( known_segments_1 )
{
	using namespace seqan;
	typedef String< Dna5 > string_t;
	typedef Iterator< const String< Dna5 > >::Type it_t;
	typedef boost::tuple< it_t, it_t > seg_t;
	typedef vector< seg_t > segment_vector;
	string_t seq = "NNNNNN";
	segment_vector known_segments;
	segment_vector unknown_segments;
	segmentize_by_known_bases(
		seq,
		make_segment_appender( known_segments ),
		make_segment_appender( unknown_segments )
	);

	BOOST_CHECK_EQUAL( int( unknown_segments.size() ), 1 );
	BOOST_CHECK( unknown_segments[0] == seg_t( seqan::begin( seq ), seqan::end( seq ) ) );
	BOOST_CHECK_EQUAL( int( known_segments.size() ), 0 );
}

BOOST_AUTO_TEST_CASE( known_segments_2 )
{
	using namespace seqan;
	typedef String< Dna5 > string_t;
	typedef Iterator< const String< Dna5 > >::Type it_t;
	typedef boost::tuple< it_t, it_t > seg_t;
	typedef vector< seg_t > segment_vector;
	string_t seq = "ACGGGT";
	segment_vector known_segments;
	segment_vector unknown_segments;
	segmentize_by_known_bases(
		seq,
		make_segment_appender( known_segments ),
		make_segment_appender( unknown_segments )
	);

	BOOST_CHECK_EQUAL( int( known_segments.size() ), 1 );
	BOOST_CHECK( known_segments[0] == seg_t( seqan::begin( seq ), seqan::end( seq ) ) );
	BOOST_CHECK_EQUAL( int( unknown_segments.size() ), 0 );
}

BOOST_AUTO_TEST_CASE( known_segments_3 )
{
	using namespace seqan;
	typedef String< Dna5 > string_t;
	typedef Iterator< const String< Dna5 > >::Type it_t;
	typedef boost::tuple< it_t, it_t > seg_t;
	typedef vector< seg_t > segment_vector;
	string_t seq = "";
	segment_vector known_segments;
	segment_vector unknown_segments;
	segmentize_by_known_bases(
		seq,
		make_segment_appender( known_segments ),
		make_segment_appender( unknown_segments )
	);

	BOOST_CHECK_EQUAL( int( known_segments.size() ), 0 );
	BOOST_CHECK_EQUAL( int( unknown_segments.size() ), 0 );
}

BOOST_AUTO_TEST_CASE( known_segments_4 )
{
	using namespace seqan;
	typedef String< Dna5 > string_t;
	typedef Iterator< const String< Dna5 > >::Type it_t;
	typedef boost::tuple< it_t, it_t > seg_t;
	typedef vector< seg_t > segment_vector;
	string_t seq = "NANCNGNTN";
	segment_vector known_segments;
	segment_vector unknown_segments;
	segmentize_by_known_bases(
		seq,
		make_segment_appender( known_segments ),
		make_segment_appender( unknown_segments )
	);

	BOOST_CHECK_EQUAL( int( known_segments.size() ), 4 );
	BOOST_CHECK( known_segments[0] == seg_t( seqan::begin( seq ) + 1, seqan::begin( seq ) + 2 ) );
	BOOST_CHECK( known_segments[1] == seg_t( seqan::begin( seq ) + 3, seqan::begin( seq ) + 4 ) );
	BOOST_CHECK( known_segments[2] == seg_t( seqan::begin( seq ) + 5, seqan::begin( seq ) + 6 ) );
	BOOST_CHECK( known_segments[3] == seg_t( seqan::begin( seq ) + 7, seqan::begin( seq ) + 8 ) );
	BOOST_CHECK_EQUAL( int( unknown_segments.size() ), 5 );
	BOOST_CHECK( unknown_segments[0] == seg_t( seqan::begin( seq ) + 0, seqan::begin( seq ) + 1 ) );
	BOOST_CHECK( unknown_segments[1] == seg_t( seqan::begin( seq ) + 2, seqan::begin( seq ) + 3 ) );
	BOOST_CHECK( unknown_segments[2] == seg_t( seqan::begin( seq ) + 4, seqan::begin( seq ) + 5 ) );
	BOOST_CHECK( unknown_segments[3] == seg_t( seqan::begin( seq ) + 6, seqan::begin( seq ) + 7 ) );
	BOOST_CHECK( unknown_segments[4] == seg_t( seqan::begin( seq ) + 8, seqan::begin( seq ) + 9 ) );
}





