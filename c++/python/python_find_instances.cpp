/** Copyright John Reid 2011
 *
 * \file Exposes find instances part of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/find_instances.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;


typedef find_instances< STEME_INDEX_MODULE_SPEC > find_instances_t;
typedef data< STEME_INDEX_MODULE_SPEC > data_t;
typedef model< STEME_INDEX_MODULE_SPEC > model_t;



void
expose_find_instances() {


    py::class_<
        find_instances_t,
        find_instances_t::ptr,
        boost::noncopyable
    > find_instances_class(
        "FindInstances",
        "Finds the instances of the model in the data (above a given Z threshold).",
        py::init<
            data_t &,
            model_t &,
            double
        >(
            (
                py::args( "data" ),
                py::args( "model" ),
                py::args( "Z_threshold" )
            ),
            "Constructor."
        )
    );
    find_instances_class.def( "__call__", &find_instances_t::descend_tree, "Find the best W-mers." );
    find_instances_class.def_readonly( "instances", &find_instances_t::instances, "Instances." );
    expose_tree_descender< find_instances_t >( find_instances_class );
}


