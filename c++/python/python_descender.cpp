/** Copyright John Reid 2011
 *
 * \file Exposes model part of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/descender.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;



void
expose_descender() {
	//
	// Stem::stats
	//
	py::class_< efficiency_stats > stats_class(
		"Stats",
		"Statistics for the STEM algorithm.",
		py::no_init
	);
	stats_class.def( "reset", &efficiency_stats::reset, "Reset to 0 counts." );
	stats_class.def_readonly( "counts", &efficiency_stats::counts, "Counts of nodes/occurrences evaluated/discarded." );

	{
		py::scope scope = stats_class;

		py::class_< efficiency_stats::count > count_class(
			"Count",
			"Counts for the STEM algorithm.",
			py::no_init
		);
		count_class.def_readonly( "evaluated", &efficiency_stats::count::evaluated, "# evaluated." );
		count_class.def_readonly( "discarded", &efficiency_stats::count::discarded, "# discarded." );
		count_class.def( "reset", &efficiency_stats::count::reset, "Reset to 0 counts." );
		count_class.add_property( "total", &efficiency_stats::count::total, "total evaluated & discarded." );
		count_class.add_property( "fraction_evaluated", &efficiency_stats::count::fraction_evaluated, "fraction evaluated out of total." );

		py::class_< efficiency_stats::node_counts > node_count_class(
			"NodeCount",
			"Count for a node in the STEM algorithm.",
			py::no_init
		);
		node_count_class.def_readonly( "node", &efficiency_stats::node_counts::node_count, "Count for nodes." );
		node_count_class.def_readonly( "occurrence", &efficiency_stats::node_counts::occurrence_count, "Counts for occurrences." );
		node_count_class.def( "reset", &efficiency_stats::node_counts::reset, "Reset to 0 counts." );

		py::class_<
			efficiency_stats::node_counts::vec,
			boost::shared_ptr< efficiency_stats::node_counts::vec >
		> node_count_vec_class( "CountVec" );
		node_count_vec_class.def( py::indexing::container_suite< efficiency_stats::node_counts::vec >() );
	}
}

