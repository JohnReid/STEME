/** Copyright John Reid 2011
 *
 * \file Exposes Data in STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/data.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;

typedef data< STEME_INDEX_MODULE_SPEC > data_t;


/**
 * Get the subsequence of given sequence with given position and length.
 */
std::string
get_data_subsequence( data_t & _data, size_t sequence, size_t position, size_t length ) {
	using namespace seqan;
	seqan_types_t::text_t & text = _data.get_text();
	const seqan_types_t::id_t _id = positionToId( text, sequence );
	return std_string_from_seqan( infix( valueById( text, _id ), position, position+length ) );
}



/** Create a new data object. */
data_t::ptr
create_data(
	seqan_types_t::index_t & index,
	bool using_dense = false,
	size_t max_W = 0
) {
	return data_t::ptr( new data_t( index, using_dense, max_W ) );
}



void
expose_data() {


	//
	// Base-resolution position-specific priors
	//
	py::class_<
		base_resolution_psp,
		boost::noncopyable,
		boost::shared_ptr< base_resolution_psp >
	> base_resolution_psp_class(
		"BaseResolutionPSP",
		"Base-resolution position-specific priors.",
		py::no_init
	);
	base_resolution_psp_class.def(
	    "adjust_for_binding_site",
	    &base_resolution_psp::adjust_for_binding_site,
        (
            py::arg( "seq" ),
            py::arg( "seq_pos" ),
            py::arg( "W" ),
            py::arg( "Z" )
        ),
	    "Adjust the PSP to account for a binding site at the given location with the given expected Z." );



	//
	// Data
	//
	py::class_<
		data_t,
		boost::noncopyable,
		data_t::ptr
	> data_class(
		"Data",
		"The data.",
        py::init<
            seqan_types_t::index_t &,
            bool,
            size_t
        >(
            (
                py::arg( "index" ),
                py::arg( "using_dense" ) = false,
                py::arg( "max_W" ) = 0
            ),
            "Constructor."
        )[
            py::with_custodian_and_ward< 1, 2 >() // the index should stay alive as long as data
        ]
	);
	data_class.add_property( "index", py::make_function( &get_index< data_t >, py::return_internal_reference<>() ), "The index." );
    data_class.def( "update_priors", &data_t::update_priors, "Update the width-specific priors after the base-resolution priors have changed." );
    data_class.def( "get_W_mer", &data_t::get_W_mer, "Get the W-mer at the given global position." );
	//data_class.def( "get_node_count", &data_t::get_node_count, "How many nodes at the given depth?" );
	data_class.def( "pos_localise", &data_t::pos_localise, "Convert a global position to a local position." );
	data_class.def( "pos_globalise", &data_t::pos_globalise, "Convert a local position to a global position." );
    data_class.def( "num_W_mers", &data_t::num_W_mers, "Number of W-mers in the data." );
    data_class.def( "num_occurrences", &data_t::get_occurrence_count, "Number of known W-mers in the data." );
	data_class.def( "seq_length", &data_t::get_sequence_length, "Get the length of given sequence." );
	data_class.def( "subsequence", get_data_subsequence, "Get the subsequence of given sequence with given position and length." );
	data_class.def_readonly( "N", &data_t::N, "Size of the data." );
	data_class.def_readonly( "base_prior", &data_t::_base_prior, "Base-resolution position-specific priors." );
	data_class.add_property( "num_sequences", &data_t::num_sequences, "Number of sequences." );
	myrrh::python::register_tuple_converter< seqan_types_t::local_pos_t >();
}


