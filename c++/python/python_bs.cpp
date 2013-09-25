/** Copyright John Reid 2011
 *
 * \file Exposes binding model parts of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/binding_model.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;







PssmStorage::matrix & get_log_prob_values( PssmStorage & log_probs ) { return log_probs.m; }


void
seed_bs_model( PssmBindingSiteModel & s, const std::string & seed, bool use_pseudo_count = true ) {
	s.seed( seqan_types_t::string_t( seed ), use_pseudo_count );
}

/// The information content of the column
double
column_IC( const Pssm & pssm, zero_order_frequencies freqs, size_t col ) {
	return calculate_column_ic( pssm.log_probs.m[ col ], freqs.dist_logs.data() );
}



void
expose_bs() {

	//
	// PssmStorage
	//
	py::class_< PssmStorage, PssmStorage::ptr > (
		"PssmStorage",
		"Holds a PSSM's probabilities, weights or log probabilities as usage requires.",
		py::init< size_t, size_t >( py::args( "W", "alphabet_size" ), "Construct a PssmStorage object." )
    )
	.def( "W", &PssmStorage::W )
	.def( "alphabet_size", &PssmStorage::alphabet_size )
	.def( "values", get_log_prob_values, "Get the values stored in the PSSM.", py::return_value_policy< py::return_by_value >() )
	//.add_property( "values", py::make_function( get_log_prob_values_py, py::with_custodian_and_ward< 0, 1 >() ), "The values stored in the PSSM." ) // this is what we would like to have
    ;




	//
	// Pssm
	//
	{
		py::class_<
			Pssm,
			Pssm::ptr,
			boost::noncopyable
		> pssm_class(
			"Pssm",
			"A PSSM.",
			py::init< PssmStorage >(
				py::arg( "log_probs" ),
				"Construct a Pssm object."
			)
		);
		py::scope scope = pssm_class;

		pssm_class.add_property( "log_probs", py::make_getter( &Pssm::log_probs ), "The PSSM's log probabilities." );
		pssm_class.def_readwrite( "num_samples", &Pssm::num_samples, "Number of samples this model was generated from." );
		pssm_class.def( "columnIC", column_IC, "The information content of the column." );
	}






	//
	// PssmBindingSiteModel
	//
	{
		py::class_<
			PssmBindingSiteModel,
			PssmBindingSiteModel::ptr,
			boost::noncopyable
		> pssm_binding_site_model_class(
			"PssmBindingSiteModel",
			"A binding site model for PSSMs.",
			py::init< PssmStorage &, double >(
				( py::arg( "pssm_log_probs" ), py::arg( "seed_pseudo_counts" )=.5 ),
				"Construct a PssmBindingSiteModel object."
			)
		);
		py::scope scope = pssm_binding_site_model_class;

		pssm_binding_site_model_class.add_property( "pssm", py::make_getter( &PssmBindingSiteModel::pssm ), "The PSSM for this binding site model." );
        pssm_binding_site_model_class.def( "seed", seed_bs_model, "Seed the model using a w-mer." );
        pssm_binding_site_model_class.def( "recalculate", &PssmBindingSiteModel::recalculate, "Recalculate the best prefixes and suffixes after the log probabilites have changed." );
		pssm_binding_site_model_class.def_readwrite( "seed_pseudo_counts", &PssmBindingSiteModel::seed_pseudo_counts, "Pseudo-counts used when seeding model." );

	}



}
