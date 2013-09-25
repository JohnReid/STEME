/*
 * steme_python_defs.h
 *
 *  Created on: 13 Aug 2011
 *      Author: john
 */

#ifndef STEME_JR_13AUG2011_STEME_PYTHON_DEFS_H_
#define STEME_JR_13AUG2011_STEME_PYTHON_DEFS_H_


#include "steme_python_pch.h"
#include <boost/python.hpp>
#include "steme_defs_pch.h"


#define STEME_INDEX_DEFAULT_TYPE 0
#define STEME_INDEX_GENOME_TYPE 1


// Choose default specification if not defined
#ifndef STEME_INDEX_MODULE_TYPE
# define STEME_INDEX_MODULE_TYPE STEME_INDEX_DEFAULT_TYPE
#endif //STEME_INDEX_MODULE_TYPE


#if STEME_INDEX_MODULE_TYPE == STEME_INDEX_DEFAULT_TYPE
# define STEME_INDEX_MODULE_SPEC ::steme::default_spec
# include <steme/spec_default.h>

#elif STEME_INDEX_MODULE_TYPE == STEME_INDEX_GENOME_TYPE
# define STEME_INDEX_MODULE_SPEC ::steme::genome_spec
# include <steme/spec_genome.h>

#else
# error Unknown index module type

#endif //STEME_INDEX_MODULE_SPEC


#include <boost/python/make_function.hpp>

#include <myrrh/python/tuple.h>
#include <myrrh/python/multi_array_to_numpy.h>
#include <myrrh/python/boost_function.h>

#include <indexing_suite/vector.hpp>
#include <indexing_suite/map.hpp>







// choose maximum Markov model order if not decided already.
// N.B. This is an upper limit so setting this to 5 will make order-4 the largest usable order
#ifndef STEME_MAX_MARKOV_MODEL_ORDER
# define STEME_MAX_MARKOV_MODEL_ORDER 7
#endif // STEME_MAX_MARKOV_MODEL_ORDER




namespace steme {

//forward decls
template< typename Spec > struct model;
template< typename Spec > struct data;


namespace python {


typedef steme_seqan_types< STEME_INDEX_MODULE_SPEC > seqan_types_t;


/// The model.
template< typename Exposed >
model< STEME_INDEX_MODULE_SPEC > &
get_model( Exposed & e ) {
	return e._model;
}


/// The data.
template< typename Exposed >
data< STEME_INDEX_MODULE_SPEC > &
get_data( Exposed & e ) {
	return e._data;
}

/// The index
template< typename Exposed >
seqan_types_t::index_t &
get_index( Exposed & e ) {
	return e.index;
}


template< typename Exposed, typename Class >
void
expose_tree_descender( Class & class_ ) {
	namespace py = boost::python;

	class_.add_property( "model", py::make_function( &get_model< Exposed >, py::return_internal_reference<>() ), "The model." );
	class_.add_property( "data", py::make_function( &get_data< Exposed >, py::return_internal_reference<>() ), "The data." );
	class_.def_readonly( "W", &Exposed::W, "The width of the model." );
	class_.def_readonly( "efficiency_statistics", &Exposed::stats, "Statistics for how many nodes evaluated." );
}





template< typename T >
bool optional_has_value( boost::optional< T > optional ) {
	return optional;
}

template< typename T >
T
optional_value( boost::optional< T > optional ) {
	if( ! optional ) {
		throw std::logic_error( "Optional has no value." );
	}
	return *optional;
}

template< typename T >
void
expose_optional( const char * optional_class_name ) {
	namespace py = boost::python;

	typedef boost::optional< T >           optional_t;

	py::class_<
		optional_t
	> optional_class(
	    optional_class_name,
        py::init< T const & >( py::arg( "value" ), "Constructor." )
    );
    optional_class.def( py::init<>() );
	optional_class.def( "__nonzero__", optional_has_value< T >, "True iff this optional has a value." );
	optional_class.add_property( "value", optional_value< T >, "The value of this optional." );
}


template< typename T >
int
cmp( T lhs, T rhs )
{
	if( lhs < rhs ) return -1;
	if( rhs < lhs ) return 1;
	return 0;
}





template< typename T >
void
expose_optional_pair( const char * optional_class_name, const char * pair_class_name ) {
	namespace py = boost::python;

	typedef optional_pair< T >                           pair_t;

	// Expose the optional type.
	expose_optional< T >( optional_class_name );

	// Expose the optional pair type.
	py::class_<
		pair_t
	> optional_pair_class( pair_class_name, py::no_init );
	optional_pair_class.def_readonly( "first", &pair_t::first, "The first member of the optional pair. Normally the value for the positive strand." );
	optional_pair_class.def_readonly( "second", &pair_t::second, "The second member of the optional pair. Normally the value for the negative strand." );
}





} // namespace steme
} // namespace python




#endif /* STEME_JR_13AUG2011_STEME_PYTHON_DEFS_H_ */
