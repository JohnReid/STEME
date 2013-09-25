/** Copyright John Reid 2011
 *
 * \file Exposes LLR p-value calculators to python.
 */

#include "steme_python_defs.h"

#include <pvalues/calculator.h>

#include <myrrh/python/boost_range.h>

#include <boost/scoped_array.hpp>

namespace py = boost::python;
using namespace pvalues;


/**
 * Create a p-value calculator that uses the FAST code.
 */
llr_pvalue_calculator::shared_ptr
py_create_fast_pvalue_calculator( py::object py_bg, size_t N, size_t Q ) {
    const size_t K = py::len( py_bg );
    boost::scoped_array< double > bg( new double[ K ] );
    std::copy(
        boost::make_transform_iterator( myrrh::python::py_seq_begin( py_bg ), myrrh::python::extract_fn< double >() ),
        boost::make_transform_iterator( myrrh::python::py_seq_end( py_bg ), myrrh::python::extract_fn< double >() ),
        bg.get()
    );
    return
        create_fast_pvalue_calculator(
            bg.get(),
            K,
            N,
            Q
        );
}



/**
 * Create a p-value calculator based on the shifted Hirji algorithm.
 */
llr_pvalue_calculator::shared_ptr
py_create_shifted_hirji_pvalue_calculator( py::object py_bg, size_t N, size_t Q ) {
    const size_t K = py::len( py_bg );
    boost::scoped_array< double > bg( new double[ K ] );
    std::copy(
        boost::make_transform_iterator( myrrh::python::py_seq_begin( py_bg ), myrrh::python::extract_fn< double >() ),
        boost::make_transform_iterator( myrrh::python::py_seq_end( py_bg ), myrrh::python::extract_fn< double >() ),
        bg.get()
    );
    return
        create_shifted_hirji_pvalue_calculator(
            bg.get(),
            K,
            N,
            Q
        );
}


/**
 * Create a p-value calculator that uses the Bejerano algorithm.
 */
boost::shared_ptr< llr_pvalue_calculator >
py_create_bejerano_pvalue_calculator( py::object py_bg, size_t max_N ) {
    const size_t K = py::len( py_bg );
    boost::scoped_array< double > bg( new double[ K ] );
    std::copy(
        boost::make_transform_iterator( myrrh::python::py_seq_begin( py_bg ), myrrh::python::extract_fn< double >() ),
        boost::make_transform_iterator( myrrh::python::py_seq_end( py_bg ), myrrh::python::extract_fn< double >() ),
        bg.get()
    );
    return
        create_bejerano_pvalue_calculator(
            bg.get(),
            K,
            max_N
        );
}



void
expose_llr_pvalues() {

    //
    // Calculate LLR's p-value
    //
    py::class_<
        llr_pvalue_calculator,
        llr_pvalue_calculator::shared_ptr,
        boost::noncopyable
    > llr_pvalue_calculator_class(
        "LLRpValueCalculator",
        "Calculates the p-values for log-likelihood ratios.",
        py::no_init
    );
    llr_pvalue_calculator_class.def( "__call__", &llr_pvalue_calculator::operator(), "Calculate the p-value of the LLR." );
    llr_pvalue_calculator_class.add_property( "maxN", &llr_pvalue_calculator::get_max_N, "The largest N this calculator can handle." );

    py::def(
        "create_bejerano_pvalue_calculator",
        py_create_bejerano_pvalue_calculator,
        (
            py::arg( "background" ),
            py::arg( "maxN" )
        ),
        "Create a LLR p-value calculator that uses Bejerano's convex algorithm."
    );

    py::def(
        "create_shifted_hirji_pvalue_calculator",
        py_create_shifted_hirji_pvalue_calculator,
        (
            py::arg( "background" ),
            py::arg( "N" ),
            py::arg( "Q" ) = 0
        ),
        "Create a LLR p-value calculator that uses our implementation of the shifted Hirji algorithm by Nagarajan."
    );

    py::def(
        "create_fast_pvalue_calculator",
        py_create_fast_pvalue_calculator,
        (
            py::arg( "background" ),
            py::arg( "N" ),
            py::arg( "Q" )
        ),
        "Create a LLR p-value calculator that uses Nagarajan's implementation of the shifted Hirji algorithm."
    );
}

