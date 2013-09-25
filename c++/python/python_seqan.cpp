/** Copyright John Reid 2011, 2013
 *
 * \file Exposes seqan parts of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/background_model.h>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <myrrh/python/boost_range.h>


namespace py = boost::python;
using namespace steme;
using namespace steme::python;

typedef std::vector< std::string > string_vector;
typedef seqan::Value< seqan::Fibre< seqan_types_t::index_t, seqan::FibreSA >::Type >::Type fibre_sa_t;
typedef seqan::Value< fibre_sa_t, 1 >::Type seq_index_t;
typedef seqan::Value< fibre_sa_t, 2 >::Type seq_offset_t;
typedef seqan::Value< seqan::Fibre< seqan_types_t::index_t, seqan::FibreChildtab >::Type >::Type text_offset_t;

typedef unsigned long long max_representing_t;
static const max_representing_t index_max_seqs = ( max_representing_t )( std::numeric_limits< seqan_types_t::seq_index_t >::max() ) + 1;
static const max_representing_t index_max_seq_length = ( max_representing_t )( std::numeric_limits< seqan_types_t::seq_offset_t >::max() ) + 1;
static const max_representing_t index_max_text_length = ( max_representing_t )( std::numeric_limits< seqan_types_t::global_pos_t >::max() ) + 1;


/**
 * Visit a node in the suffix tree (array).
 */
void
visit_node( seqan_types_t::top_down_it it, seqan::Size< seqan_types_t::index_t >::Type max_depth ) {
    using namespace seqan;
    const bool reached_max_depth = 0 != max_depth && repLength( it ) >= max_depth;
    if( ! reached_max_depth && goDown( it ) ) {
        visit_node( it, max_depth );
        while( goRight( it ) ) {
            visit_node( it, max_depth );
        }
    }
}

/**
 * Visit all the nodes in a tree. This can be useful to build the suffix tree ahead of doing any timings.
 */
void
visit_tree( seqan_types_t::index_t & index, seqan::Size< seqan_types_t::index_t >::Type max_depth ) {
    visit_node( seqan_types_t::top_down_it( index ), max_depth );
}



/**
 * Read a FASTA file into a string set.
 */
size_t
read_fasta_( const char * fasta_filename, seqan_types_t::string_set_t & sequences, string_vector & ids )
{
    using namespace seqan;

    ::std::ifstream f( fasta_filename, ::std::ios_base::in | ::std::ios_base::binary );
    size_t num_bases = 0;
    if( ! f ) {
        throw std::logic_error( MYRRH_MAKE_STRING( "Could not open: \""<<fasta_filename<<"\"" ) );
    } else {
        //scope_progress_timer timer( std::cout, MYRRH_MAKE_STRING( "Reading fasta file \"" << fasta_filename << "\"" ) );
        String< char > meta;
        seqan_types_t::string_t str;
        while( ! _streamEOF( f ) ) {
            readMeta( f, meta, Fasta() );
            read( f, str, Fasta() );
            if( length( str ) == 0 && length( meta ) == 0 ) continue;
            num_bases += length( str );
            std::string tmp_meta;
            assign( tmp_meta, meta );
            ids.push_back( tmp_meta );
            appendValue( sequences, str );
        }
        //std::cout << "Read " << length( sequences ) << " sequences with a total of " << num_bases << " bases\n";
        f.close();
    }
    return num_bases;
}

/** Read a FASTA file into a seqan string set. */
py::tuple
read_fasta( const char * fasta_filename ) {
    boost::shared_ptr< string_vector > ids( new string_vector );
    seqan_types_t::string_set_ptr_t sequences( new seqan_types_t::string_set_t );
    const size_t num_bases = read_fasta_( fasta_filename, *sequences, *ids );
    return py::make_tuple( num_bases, sequences, ids );
}

size_t
string_set_len( seqan_types_t::string_set_ptr_t string_set ) {
    return seqan::length( *string_set );
}

size_t
string_len( seqan_types_t::string_set_ptr_t string_set, size_t string_index ) {
    using namespace seqan;
    return length( getValueById( *string_set, string_index ) );
}


void
append_string( seqan_types_t::string_set_ptr_t string_set, std::string s ) {
    using namespace seqan;
    seqan_types_t::string_t seqan_str;
    assign( seqan_str, s );
    appendValue( *string_set, seqan_str );
}

/** Build an index from a string set. */
seqan_types_t::index_ptr_t
build_index( const seqan_types_t::string_set_t & sequences )
{
    if( index_max_seqs < length( sequences ) ) {
        throw std::length_error( "Too many sequences for this index type." );
    }
    max_representing_t total = 0;
    BOOST_FOREACH( const seqan_types_t::string_t & seq, sequences ) {
        if( index_max_seq_length < length( seq ) ) {
            throw std::length_error( "Sequence too long for this index type." );
        }
        total += length( seq );
    }
    if( total > index_max_text_length ) {
        throw std::length_error( "Total sequence length too big for this index type." );
    }
    return seqan_types_t::index_ptr_t( new seqan_types_t::index_t( sequences ) );
}

/** Save an index to a file. */
void
save_index( seqan_types_t::index_ptr_t index, const std::string & filename )
{
    if( ! seqan::save( *index, filename.c_str() ) ) {
        throw std::logic_error(
            MYRRH_MAKE_STRING( "Could not save index to " << filename ) );
    }
}

/** Load an index from a file. */
seqan_types_t::index_ptr_t
load_index( const std::string & filename )
{
    seqan_types_t::index_ptr_t result( new seqan_types_t::index_t );
    if( ! seqan::open( *result, filename.c_str() ) ) {
        throw std::logic_error(
            MYRRH_MAKE_STRING( "Could not load index from " << filename ) );
    }
    return result;
}

const seqan_types_t::string_set_t &
get_index_text( seqan_types_t::index_t & index ) {
    return seqan::indexText( index );
}

seqan_types_t::string_set_limits_ptr
create_string_set_limits( seqan_types_t::string_set_ptr_t string_set ) {
    seqan_types_t::string_set_limits_ptr result( new seqan_types_t::string_set_limits_t );
    *result = seqan::stringSetLimits( *string_set );
    return result;
}

std::string
get_representative( seqan_types_t::top_down_it it ) {
    std::string result;
    seqan::append( result, seqan::representative( it ) );
    return result;
}

seqan::Size< seqan_types_t::index_t >::Type
count_occurrences( seqan_types_t::top_down_it it ) {
    return seqan::countOccurrences( it );
}

bool
go_down( seqan_types_t::top_down_it & it ) {
    return seqan::goDown( it );
}

bool
go_down_base( seqan_types_t::top_down_it & it, char c ) {
    return seqan::goDown( it, seqan_types_t::alphabet_t( c ) );
}

bool
go_right( seqan_types_t::top_down_it & it ) {
    return seqan::goRight( it );
}



template< typename SeqanStr >
struct seqan_to_py_string
{
    static
    PyObject *
    convert( SeqanStr const & seq_str )
    {
        const std::string py_str = std_string_from_seqan( seq_str );
        return PyString_FromString( py_str.c_str() );
    }
};

template< typename SeqanStr >
void
register_seqan_string_to_python() {
    py::to_python_converter<
        SeqanStr,
        seqan_to_py_string< SeqanStr >
    >();
}


template< typename T >
void dummy() {
    T t = "";
}


struct string_pickle_suite : boost::python::pickle_suite
{
    static
    boost::python::tuple
    getinitargs( const seqan_types_t::string_t & x )
    {
        return boost::python::tuple();
    }

    static
    boost::python::tuple
    getstate( const seqan_types_t::string_t & x )
    {
        std::string s;
        assign( s, x );
        return boost::python::make_tuple( s );
    }

    static
    void
    setstate( seqan_types_t::string_t & x, boost::python::tuple state )
    {
        std::string s = boost::python::extract< std::string >( state[ 0 ] );
        assign( x, s );
    }
};


struct string_set_pickle_suite : boost::python::pickle_suite
{
    static
    boost::python::tuple
    getinitargs( const seqan_types_t::string_set_t & x )
    {
        return boost::python::make_tuple();
    }

    static
    boost::python::tuple
    getstate( const seqan_types_t::string_set_t & x )
    {
        boost::python::list sequence;
        BOOST_FOREACH( const seqan_types_t::string_t & s, x ) {
            sequence.append( s );
        }
        return boost::python::tuple( sequence );
    }

    static
    void
    setstate( seqan_types_t::string_set_t & x, boost::python::tuple state )
    {
        BOOST_FOREACH(
            const seqan_types_t::string_t & s,
            myrrh::python::make_boost_range< seqan_types_t::string_t >( state )
        ) {
            seqan::appendValue( x, s );
        }
    }
};


struct string_vec_pickle_suite : boost::python::pickle_suite
{
    static
    boost::python::tuple
    getinitargs( const string_vector & x )
    {
        return boost::python::make_tuple();
    }

    static
    boost::python::tuple
    getstate( const string_vector & x )
    {
        boost::python::list sequence;
        BOOST_FOREACH( const std::string & s, x ) {
            sequence.append( s );
        }
        return boost::python::tuple( sequence );
    }

    static
    void
    setstate( string_vector & x, boost::python::tuple state )
    {
        BOOST_FOREACH(
            const std::string & s,
            myrrh::python::make_boost_range< std::string >( state )
        ) {
            x.push_back( s );
        }
    }
};


void
expose_seqan() {

    py::scope().attr( "index_type" ) = STEME_INDEX_DEFAULT_TYPE;

    register_seqan_string_to_python< seqan_types_t::infix_t >();

    py::scope().attr( "index_max_seqs" ) = index_max_seqs;
    py::scope().attr( "index_max_seq_length" ) = index_max_seq_length;
    py::scope().attr( "index_max_text_length" ) = index_max_text_length;


    py::class_<
        string_vector,
        boost::shared_ptr< string_vector >,
        boost::noncopyable
    > string_vec_class( "StringVec", "A sequence of strings" );
    string_vec_class.def( py::indexing::container_suite< string_vector >() );
    string_vec_class.def_pickle( string_vec_pickle_suite() );

    py::def( "read_fasta", read_fasta, "Read a FASTA file into a seqan string set." );
    py::def( "build_index", build_index, "Build a seqan index from a string set." );
    py::def( "save_index", save_index, "Save a seqan index to a file." );
    py::def( "load_index", load_index, "Load a seqan index from a file." );

    //
    // StringSet
    //
    py::class_<
        seqan_types_t::string_t,
        boost::noncopyable,
        boost::shared_ptr< seqan_types_t::string_t >
    > string_class(
        "String",
        "Wrapper for C++ seqan string."
    );
    string_class.def_pickle( string_pickle_suite() );

    //
    // StringSet
    //
    py::class_<
        seqan_types_t::string_set_t,
        boost::noncopyable,
        seqan_types_t::string_set_ptr_t
    > string_set_class(
        "StringSet",
        "Wrapper for C++ seqan string set."
    );
    string_set_class.def( "__len__", string_set_len, "The number of strings in the string set." );
    string_set_class.def( "string_length", string_len, "The length of one of the strings in the string set." );
    string_set_class.def( "append", append_string, "Append a string to the string set." );
    string_set_class.def_pickle( string_set_pickle_suite() );
    //py::register_ptr_to_python< seqan_types_t::string_set_ptr_t >();




    //
    // Index
    //
    py::class_<
        seqan_types_t::index_t,
        seqan_types_t::index_ptr_t,
        boost::noncopyable
    > index_class(
        "Index",
        "Wrapper for C++ seqan index."
    );
    index_class.def(
        "visit",
        visit_tree,
        ( py::arg( "max_depth" ) = 0 ),
        "Visit all the nodes in a tree. This can be useful to build the suffix tree ahead of doing any timings."
    );
    index_class.def( "text", get_index_text, "The text of this index.", py::return_internal_reference<>() );



    //
    // Iterator
    //
    py::class_< seqan_types_t::top_down_it > it_class(
        "Iterator",
        "Iterator into text index.",
        py::init< seqan_types_t::index_t & >( py::arg( "index" ), "Constructor." )[ py::with_custodian_and_ward< 1, 2 >() ]
    );
    it_class.add_property( "representative", get_representative, "A representative substring of this iterator." );
    it_class.add_property( "num_occurrences", count_occurrences, "Number of occurrences of the prefix this iterator represents." );
    it_class.def( "goDown", go_down, "Iterates down one edge or a path in a tree." );
    it_class.def( "goDownBase", go_down_base, "Iterates down one edge or a path in a tree which starts with given base." );
    it_class.def( "goRight", go_right, "Iterates to the next sibling in a tree." );


}
