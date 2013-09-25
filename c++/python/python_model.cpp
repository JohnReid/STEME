/** Copyright John Reid 2011
 *
 * \file Exposes model part of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/model.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;


typedef model< STEME_INDEX_MODULE_SPEC > model_t;
typedef data< STEME_INDEX_MODULE_SPEC > data_t;



namespace steme {
extern size_t max_samples_in_E_value_calc;
}



/// Get prior number of sites
double
get_prior_num_sites( model_t & model ) {
	return model._prior.num_sites;
}


/// Set prior number of sites
void
set_prior_num_sites( model_t & model, double prior_num_sites ) {
	model._prior.num_sites = prior_num_sites;
}

/// The background model
model_t::bg_model_t &
get_bg( model_t & _model ) {
    return _model.bg;
}

/// The data
data_t &
get_model_data( model_t & _model ) {
    return _model._data;
}

/// Get lambda
double
get_lambda( model_t & model ) {
	return model._prior.lambda.get_lambda();
}


/// Set lambda
void
set_lambda( model_t & model, double lambda ) {
	model._prior.lambda.set_lambda( lambda );
}



void
expose_model()
{
	//
	// Model
	//
	py::class_<
		model_t,
		model_t::ptr,
		boost::noncopyable
	> model_class(
		"Model",
		"Model of the data.",
        py::init<
            data_t &,
            model_t::bs_model_t &,
            model_t::bg_model_t &,
            double
        >(
            (
                py::arg( "data" ),
                py::arg( "bs" ),
                py::arg( "bg" ),
                py::arg( "_lambda" ) = 0.
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
          >() // the data, binding site model and background model should stay alive as long as the model
        ]
	);
    model_class.def( "copy", &model_t::copy, "Copy this model." );
    model_class.add_property( "bs", py::make_getter( &model_t::bs ), "The binding site model." );
    model_class.add_property( "bg", py::make_function( &get_bg, py::return_internal_reference<>() ), "The background model." );
    model_class.add_property( "data", py::make_function( &get_model_data, py::return_internal_reference<>() ), "The data." );
	model_class.add_property( "W", &model_t::W, "Width of the motif." );
    model_class.add_property( "num_W_mers", &model_t::num_W_mers, "Number of W-mers in the data." );
    model_class.add_property( "num_occurrences", &model_t::num_occurrences, "Number of known W-mers in the data." );
	model_class.add_property( "prior_num_sites", get_prior_num_sites, set_prior_num_sites, "Prior estimate of number of sites." );
    model_class.add_property( "lambda_", get_lambda, set_lambda, "lambda, the prior probability of a W-mer being a binding site." );
    model_class.def( "set_lambda_for_sites", &model_t::set_lambda_for_sites, "Set lambda for the given number of sites." );
    model_class.def(
        "calculate_Z",
        &model_t::calculate_Z,
        (
            py::arg( "global_pos" ),
            py::arg( "rev_comp" ),
            py::arg( "g" ) = boost::optional< double >()
        ),
        "Calculate Z."
    );
}
