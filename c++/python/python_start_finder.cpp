/** Copyright John Reid 2011
 *
 * \file Exposes model part of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/start_finder.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;


typedef start_finder< DEFAULT_FIND_WMER_STORAGE, STEME_INDEX_MODULE_SPEC > start_finder_t;
typedef data< STEME_INDEX_MODULE_SPEC > data_t;
typedef model< STEME_INDEX_MODULE_SPEC > model_t;
typedef significance< STEME_INDEX_MODULE_SPEC > significance_t;


void
find_starts_for_seed( start_finder_t & start_finder, const std::string & seed ) {
	start_finder.find_starts_for_seed( seqan_types_t::string_t( seed ) );
}


/// Create a start counter
start_finder_t::start_counter_for_width
create_start_counter(
	data_t & data,
	size_t W
) {
	return start_finder_t::start_counter_for_width( data, W );
}



void
expose_start_finder() {

    py::def(
        "edit_distance",
        edit_distance< std::string, std::string >,
        (
            py::arg( "string1" ),
            py::arg( "string2" )
        ),
        "The edit distance between two strings."
    );


    py::def(
        "shifted_edit_distance",
        shifted_edit_distance< std::string, std::string >,
        (
            py::arg( "string1" ),
            py::arg( "string2" )
        ),
        "The edit distance between two strings allowing a shift as one edit."
    );


    py::def(
        "rev_comp_shifted_edit_distance",
        rev_comp_shifted_edit_distance< std::string, std::string >,
        (
            py::arg( "string1" ),
            py::arg( "string2" )
        ),
        "The edit distance between two string allowing reverse complements and a shift as one edit."
    );


    py::def(
        "create_start_counter",
        create_start_counter,
        (
            py::arg( "data" ),
            py::arg( "W" )
        ),
        py::with_custodian_and_ward_postcall< 0, 1 >() // the data should stay alive as long as counter
    );


	//
	// Find starting points
	//
	py::class_<
		start_finder_t::start,
		boost::shared_ptr< start_finder_t::start >
	> start_class(
	    "Start",
	    "A starting point for EM.",
	    py::no_init
	);
	start_class.def_readwrite( "score", &start_finder_t::start::score, "The score's start." );
	start_class.def_readwrite( "num_sites", &start_finder_t::start::num_sites, "The number of sites (W-mers) aligned to make the start." );
	start_class.def_readwrite( "seed", &start_finder_t::start::w_mer, "The W-mer that seeded this start." );
	start_class.def_readwrite( "best_w_mers", &start_finder_t::start::best_w_mers, "The best W-mers aligned to make the start." );
	start_class.def_readwrite( "model", &start_finder_t::start::_model, "The model used to evaluate the start." );

	py::class_<
		start_finder_t,
		boost::noncopyable,
		start_finder_t::ptr
	> start_finder_class(
	    "FindStarts",
        "Finds the best starting points for EM.",
        py::init<
            data_t &,
            model_t &,
            significance_t &,
            size_t,
            size_t,
            double,
            size_t,
            size_t,
            size_t
        >(
            (
                py::arg( "data" ),
                py::arg( "_model" ),
                py::arg( "_significance" ),
                py::arg( "min_num_sites" ),
                py::arg( "max_num_sites" ),
                py::arg( "candidate_starts_factor" ),
                py::arg( "num_to_find" ) = 1,
                py::arg( "speed_up" ) = 0,
                py::arg( "first_to_examine" ) = 0
            ),
            "Constructor."
        )[
          py::with_custodian_and_ward<
              1,
              2,
              py::with_custodian_and_ward<
                  1,
                  3,
                  py::with_custodian_and_ward<
                      1,
                      4
                  >
              >
          >() // the data, model and significance should stay alive as long as the start finder
        ]
	);
	start_finder_class.def( "find_starts", &start_finder_t::find_starts, "Find the best starting points." );
	start_finder_class.def( "find_starts_for_seed", find_starts_for_seed, "Find the starting points just for the given seed." );
	start_finder_class.def( "register_callback", &start_finder_t::register_callback, "Add a callback for when a start is examined." );
	start_finder_class.add_property( "model", py::make_function( &get_model< start_finder_t >, py::return_internal_reference<>() ), "The model." );
	start_finder_class.add_property( "data", py::make_function( &get_data< start_finder_t >, py::return_internal_reference<>() ), "The data." );
	start_finder_class.def_readonly( "W", &start_finder_t::W, "The width of the model." );
	start_finder_class.def_readonly( "start_counter", &start_finder_t::start_counter, "The index of the start that we are looking at." );
	start_finder_class.def_readonly( "starts_examined", &start_finder_t::starts_examined, "How many starts we have examined." );
	start_finder_class.def_readwrite( "speed_up", &start_finder_t::speed_up, "If >0, we only examine every n'th start." );
	start_finder_class.def_readwrite( "first_to_examine", &start_finder_t::first_to_examine, "The index of the first start to examine." );
	start_finder_class.def_readonly( "best_starts", &start_finder_t::best_starts, "Best starting points (partitioned)." );
	start_finder_class.def_readonly( "efficiency_statistics", &start_finder_t::stats, "Statistics for how many nodes evaluated." );
	{
		py::scope start_finder_scope = start_finder_class;
		myrrh::python::def_function<
			void(
			    start_finder_t::eval_vec_ptr,
				const std::string &,
				int,
				double
			)
		>(
			"start_examined_callback",
			"A callback function for when a start is examined."
		);

		py::class_<
			start_finder_t::start_counter_for_width
		> start_counter_class(
			"StartCounter",
			"Counts starts of a given width.",
			py::no_init
		);
		start_counter_class.def( "__call__", &start_finder_t::start_counter_for_width::operator(), "Count the starts of a given width." );
	}

    py::class_<
        start_finder_t::start_vec,
        boost::shared_ptr< start_finder_t::start_vec >
    > start_vec_class( "StartVec" );
    start_vec_class.def( py::indexing::container_suite< start_finder_t::start_vec >() );

    py::class_<
        start_finder_t::start_map,
        start_finder_t::start_map_ptr
    > start_map_class( "StartMap" );
    start_map_class.def( py::indexing::container_suite< start_finder_t::start_map >() );
}


