/** Copyright John Reid 2011
 *
 * \file Exposes significance calculator to python.
 */

#include "steme_python_defs.h"

#include <steme/significance.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;


typedef significance< STEME_INDEX_MODULE_SPEC > significance_t;


void
expose_significance() {

    //
    // Calculate model's significance
    //
    py::class_<
        significance_t,
        significance_t::ptr,
        boost::noncopyable
    > significance_class(
        "Significance",
        "Calculates the significance of different models.",
        py::init<
            significance_t::data_t &,
            zero_order_frequencies &,
            pvalues::llr_pvalue_calculator::shared_ptr
        >(
            (
                py::args( "data" ),
                py::args( "freqs" ),
                py::args( "p_value_calculator" )
            ),
            "Constructor."
        )
    );
    significance_class.def( "log_E_value", &significance_t::log_E_value, "The logarithm of the E-value of the model." );
    significance_class.def( "log_product_p_values", &significance_t::log_product_p_values, "The logarithm of the product of p-values of the model." );
}

