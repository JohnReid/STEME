/** Copyright John Reid 2011
 *
 * \file Exposes Markov model parts of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/background_model.h>
#include <steme/data.h>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;

typedef data< STEME_INDEX_MODULE_SPEC > data_t;









template< typename CompleteMarkovModel >
likelihoods_vec_vec_ptr
_calculate_likelihoods( CompleteMarkovModel & model, data_t & _data ) {
	likelihoods_vec_vec_ptr result( new likelihoods_vec_vec_t );
	calculate_likelihoods( model, _data.get_text(), *result );
	return result;
}


template< typename CompleteMarkovModel >
boost::tuple< typename CompleteMarkovModel::ptr, zero_order_frequencies >
_create_markov_model( data_t & _data, double pseudo_count ) {
	typename CompleteMarkovModel::ptr result( new CompleteMarkovModel() );
	return boost::make_tuple( result, build_model_from_text( *result, _data.get_text(), pseudo_count ) );
}

template< typename CompleteMarkovModel >
boost::tuple< typename CompleteMarkovModel::ptr, zero_order_frequencies >
_create_markov_model_from_index( seqan_types_t::index_t & index, double pseudo_count ) {
	typename CompleteMarkovModel::ptr result( new CompleteMarkovModel() );
	return boost::make_tuple( result, build_model_from_index( *result, index, pseudo_count ) );
}



template< size_t Order >
void
register_markov_model_fns() {
	typedef complete_markov_model< Order, 4, double > model_t;

	py::class_<
		model_t,
		boost::noncopyable,
		typename model_t::ptr
	> markov_model_class(
		MYRRH_MAKE_STRING( "MarkovModelOrder" << Order ).c_str(),
		MYRRH_MAKE_STRING( "A " << Order << "-order Markov model." ).c_str(),
		py::no_init
	);
	markov_model_class.def(
		"calculate_likelihoods",
		_calculate_likelihoods< model_t >,
		"Calculate the likelihoods of the sequences."
	);

	py::def(
		MYRRH_MAKE_STRING( "create_markov_model_order_" << Order ).c_str(),
		_create_markov_model< model_t >,
		"Construct a Markov model that models the sequences and also calculate the zero order frequencies."
	);

	py::def(
		MYRRH_MAKE_STRING( "create_markov_model_order_from_index_" << Order ).c_str(),
		_create_markov_model_from_index< model_t >,
		"Construct a Markov model that models the sequences and also calculate the zero order frequencies."
	);

	py::def(
		"W_mer_log_likelihood",
		W_mer_log_likelihood,
		"The likelihood of the W-mer at the given offset."
	);

	myrrh::python::register_tuple_converter< boost::tuple< typename model_t::ptr, zero_order_frequencies > >();
}


/// Registers a Markov model of given order, n
#define REGISTER_MARKOV_MODEL_OF_ORDER( z, n, text ) register_markov_model_fns< n >();

void
expose_markov() {

    BOOST_PP_REPEAT_FROM_TO( 0, STEME_MAX_MARKOV_MODEL_ORDER, REGISTER_MARKOV_MODEL_OF_ORDER, )

}

