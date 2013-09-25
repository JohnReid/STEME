/** Copyright John Reid 2011
 *
 * \file Exposes expectation-maximization part of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/em.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;


typedef EM_descender< STEME_INDEX_MODULE_SPEC > em_t;
typedef data< STEME_INDEX_MODULE_SPEC > data_t;
typedef model< STEME_INDEX_MODULE_SPEC > model_t;



/** Return an object that can do the EM algorithm. */
em_t::ptr
create_EM_descender(
	data_t & data,
	model_t & _model,
	double epsilon,
	double wnsites,
	double log_delta_factor = std::log( .9 ),
	double log_delta_adjustment_max = std::log( 5. )
) {
	return em_t::ptr( new em_t( data, _model, epsilon, wnsites ) );
}



void
expose_em() {

	py::def(
		"create_EM_descender",
		create_EM_descender,
		(
			py::arg( "data" ),
			py::arg( "model" ),
			py::arg( "epsilon" ),
			py::arg( "wnsites" ),
			py::arg( "log_delta_factor" ) = std::log( .9 ),
			py::arg( "log_delta_adjustment_max" ) = std::log( 5. )
		),
		py::with_custodian_and_ward_postcall<
			0,
			1,
			py::with_custodian_and_ward_postcall< 0, 2 >
		>() // the data and model should stay alive as long as descender
	);


	//
	// EM
	//
	py::class_<
		em_t,
		em_t::ptr
	> em_descender_class( "EM", py::no_init );
	expose_tree_descender< em_t >( em_descender_class );

	py::scope scope = em_descender_class;

	expose_optional_pair< double >( "OptionalDouble", "ZPair" );

	em_descender_class.def( "do_iteration", &em_t::do_iteration, "Perform one iteration of EM." );
	em_descender_class.def( "get_Z", &em_t::get_Z, "Evaluation (Z) at given global position." );
	em_descender_class.def_readwrite( "log_delta_factor", &em_t::log_delta_factor, "The amount by which we change log delta when we see the estimation of c_wb is too low." );
	em_descender_class.def_readwrite( "log_delta_adjustment_max", &em_t::log_delta_adjustment_max, "The largest adjustment to log delta we allow." );
	em_descender_class.def_readwrite( "log_delta_adjustment", &em_t::log_delta_adjustment, "An adjustment we make to log delta for the next iteration after seeing the actual value of c_wb. This can change every iteration." );
	em_descender_class.def_readwrite( "llr_threshold", &em_t::llr_threshold, "Threshold on log likelihood ratio." );
	em_descender_class.def_readwrite( "epsilon", &em_t::epsilon, "Allowed error in f_wb." );
	em_descender_class.def_readwrite( "wnsites", &em_t::wnsites, "Weight on number of sites. Used when updating lambda." );
	em_descender_class.def_readwrite( "using_sparse_Z", &em_t::using_sparse_Z, "Use sparse data structure to store Z." );
	em_descender_class.def_readonly( "expected_sites", &em_t::expected_sites, "Expected number of sites." );
	em_descender_class.def_readonly( "LL", &em_t::LL, "Expected log likelihood." );
	em_descender_class.def_readonly( "efficiency_statistics", &em_t::stats, "Statistics for how many nodes evaluated." );
	em_descender_class.def_readonly( "lambda_ratios", &em_t::lambda_ratios, "The lambda ratios. These tell us whether our threshold was high enough." );
	em_descender_class.def_readonly( "num_Z_non_zero", &em_t::num_Z_non_zero, "# of non-zero Z." );
	em_descender_class.def_readonly( "num_Z_large", &em_t::num_Z_large, "# of large Z." );
	em_descender_class.def_readonly( "num_Z_normalised", &em_t::num_Z_normalised, "# of Z that were normalised." );
}
