/** Copyright John Reid 2011
 *
 * \file Exposes find best W-mers part of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/find_best_w_mers.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;


typedef find_best_w_mers< DEFAULT_FIND_WMER_STORAGE, STEME_INDEX_MODULE_SPEC > find_best_w_mers_t;
typedef data< STEME_INDEX_MODULE_SPEC > data_t;
typedef model< STEME_INDEX_MODULE_SPEC > model_t;
typedef w_mer_evaluation< STEME_INDEX_MODULE_SPEC > eval_t;




/** Return an object that can find the best W-mers. */
find_best_w_mers_t::ptr
create_best_w_mer_finder(
	data_t & data,
	model_t & _model,
	size_t num_to_find
) {
	return find_best_w_mers_t::ptr( new find_best_w_mers_t( data, _model, num_to_find ) );
}


eval_t::vec_ptr
get_best_non_overlapping( eval_t::vec & instances, int W, bool already_sorted_by_position ) {
    eval_t::vec_ptr result( new eval_t::vec );
    output_non_overlapping_instances(
        instances,
        W,
        std::back_inserter( *result ),
        already_sorted_by_position
    );
    return result;
}






/** Calculate the overlap in the two vectors. */
unsigned
py_calculate_overlap( py::object evals1, size_t W1, py::object evals2, size_t W2 ) {
    using myrrh::python::make_extract_iterator;
    using myrrh::python::py_seq_begin;
    using myrrh::python::py_seq_end;

    return calculate_overlap(
        make_extract_iterator< eval_t >( py_seq_begin( evals1 ) ),
        make_extract_iterator< eval_t >( py_seq_end( evals1 ) ),
        W1,
        make_extract_iterator< eval_t >( py_seq_begin( evals2 ) ),
        make_extract_iterator< eval_t >( py_seq_end( evals2 ) ),
        W2
    );
}

template< typename FindBestWMers >
py::class_<
    FindBestWMers,
    typename FindBestWMers::ptr
>
expose_find_best_w_mers( const char * name ) {
    //
    // Find best W-mers
    //
    py::class_<
        FindBestWMers,
        typename FindBestWMers::ptr
    > find_best_w_mers_class(
        name,
        py::init<
            data_t &,
            model_t &,
            size_t
        >(
            (
                py::arg( "data" ),
                py::arg( "model" ),
                py::arg( "num_to_find" )
            ),
            "Constructor."
        )[
            py::with_custodian_and_ward<
                1,
                2,
                py::with_custodian_and_ward< 1, 3 >
            >() // the data and model should stay alive as long as best W-mer finder
        ]

    );
    find_best_w_mers_class.def( "__call__", &FindBestWMers::descend_tree, "Find the best W-mers." );
    find_best_w_mers_class.def_readonly( "best_w_mers", &FindBestWMers::best_w_mers, "Best W-mers." );
    find_best_w_mers_class.def_readonly( "has_overlapping", &FindBestWMers::has_overlapping, "Do any of the W-mers overlap?" );
    find_best_w_mers_class.def(
        "update_model",
        &FindBestWMers::update_model_with_W_mers,
        "Update the model to reflect the distribution of the W-mers.",
        ( py::arg( "num_sites" ), py::arg( "use_pseudo_counts" ) = true )
    );
    expose_tree_descender< FindBestWMers >( find_best_w_mers_class );

    return find_best_w_mers_class;
}




void
expose_find_best_w_mers() {

	py::def(
		"create_best_w_mer_finder",
		create_best_w_mer_finder,
		(
			py::arg( "data" ),
			py::arg( "model" ),
			py::arg( "num_to_find" )
		),
		py::with_custodian_and_ward_postcall<
			0,
			1,
			py::with_custodian_and_ward_postcall< 0, 2 >
		>() // the data and model should stay alive as long as best W-mer finder
	);

    py::def(
        "calculate_overlap",
        py_calculate_overlap,
        (
            py::arg( "evals1" ),
            py::arg( "W1" ),
            py::arg( "evals2" ),
            py::arg( "W2" )
        )
    );
    py::def(
        "erase_overlapping_from_score_sorted",
        erase_overlapping_from_score_sorted< STEME_INDEX_MODULE_SPEC >,
        (
            py::arg( "instances" ),
            py::arg( "W" ),
            py::arg( "max_size" ) = -1
        ),
        "Erase overlapping instances from a score sorted vector of instances."
    );
    py::def( "sort_instances_by_score", sort_instances_by_score< eval_t::vec > );
    py::def( "sort_instances_by_position", sort_instances_by_position< eval_t::vec > );

//    py::class_<
//        eval_t::vec,
//        eval_t::vec_ptr
//    > evaluations_vec_class( "EvaluationsVec" );
//    evaluations_vec_class.def( py::indexing::container_suite< eval_t::vec >() );

    expose_find_best_w_mers< find_best_w_mers< store_using_multi_index, STEME_INDEX_MODULE_SPEC > >( "FindBestWMersMultiIndex" );
    expose_find_best_w_mers< find_best_w_mers< store_in_set           , STEME_INDEX_MODULE_SPEC > >( "FindBestWMersSet" );
    expose_find_best_w_mers< find_best_w_mers< store_in_sorted_vector , STEME_INDEX_MODULE_SPEC > >( "FindBestWMersSortedVec" );

    {
        py::scope find_best_w_mers_scope = expose_find_best_w_mers< find_best_w_mers_t >( "FindBestWMers" );

        py::class_<
            eval_t
        > w_mer_evaluation_class(
            "Evaluation",
            py::init< double, eval_t::global_pos_t, bool >(
                (
                    py::arg( "Z" ),
                    py::arg( "global_pos" ),
                    py::arg( "rev_comp" )
                ),
                "Constructor."
            )
        );
        w_mer_evaluation_class.def( "__cmp__", cmp< eval_t >, "Comparison." );
        w_mer_evaluation_class.def_readonly( "Z", &eval_t::Z, "Expected value of Z for this W-mer." );
        w_mer_evaluation_class.def_readonly( "global_pos", &eval_t::global_pos, "Global position." );
        w_mer_evaluation_class.def_readonly( "rev_comp", &eval_t::rev_comp, "The reverse complement of the W-mer or not." );
    }

    py::class_<
        eval_t::vec,
        eval_t::vec_ptr
    > instance_vec_class( "InstanceVec" );
    instance_vec_class.def(
        py::indexing::container_suite< eval_t::vec >()
    );
    instance_vec_class.def( "sort_by_position", sort_instances_by_position< eval_t::vec >, "Sort the instances by position." );
    instance_vec_class.def(
        "get_best_non_overlapping",
        get_best_non_overlapping,
        (
            py::arg( "instances" ),
            py::arg( "W" ),
            py::arg( "already_sorted_by_position" ) = false
        ),
        "Get the best instances that don't overlap."
    );
    instance_vec_class.def(
        "do_instances_overlap",
        do_instances_overlap< eval_t::vec >,
        (
            py::arg( "instances" ),
            py::arg( "W" ),
            py::arg( "already_sorted_by_position" ) = false
        ),
        "Do the instances overlap?"
    );
}


