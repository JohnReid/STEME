/** Copyright John Reid 2011, 2012
 *
 * \file Exposes utility parts of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/defs.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;

#ifdef USING_GOOGLE_PROFILER
#include <google/profiler.h>
#endif //USING_GOOGLE_PROFILER




/// Useful to set break points in.
void dummy_fn() {
	int i = 0;
}




void
expose_utility() {

#ifndef NDEBUG
	py::scope().attr( "_debug" ) = true;
	std::cout << "WARNING: Debug version of _stempy module loaded. If you did not intend this then check your configuration!" << std::endl;
#else //_DEBUG
	py::scope().attr( "_debug" ) = false;
#endif //_DEBUG

    py::scope().attr( "INDEX_DEFAULT" ) = STEME_INDEX_DEFAULT_TYPE;
    py::scope().attr( "INDEX_GENOME" ) = STEME_INDEX_GENOME_TYPE;

	// expose numpy stuff
	import_array();
	myrrh::python::expose_man_fns();
	myrrh::python::expose_converters< double >();

	py::def( "_dummy_fn", dummy_fn, "An empty function that does nothing. It is useful to set breakpoints in." );


#ifdef USING_GOOGLE_PROFILER
    py::scope().attr( "_has_google_profiler" ) = true;
	py::def( "__google_profiler_start", ProfilerStart, "Start the google profiler." );
	py::def( "__google_profiler_stop", ProfilerStop, "Stop the google profiler." );
#else
    py::scope().attr( "_has_google_profiler" ) = false;
#endif //USING_GOOGLE_PROFILER


	//
	// double vector
	//
	py::class_<
		double_vec,
		boost::shared_ptr< double_vec >
	> double_vec_class( "DoubleVec" );
	double_vec_class.def( py::indexing::container_suite< double_vec >() );


	//
	// vector of double vector
	//
	py::class_<
		double_vec_vec,
		boost::noncopyable,
		boost::shared_ptr< double_vec_vec >
	> double_vec_vec_class( "DoubleVecVec" );
	double_vec_vec_class.def( py::indexing::container_suite< double_vec_vec >() );

}
