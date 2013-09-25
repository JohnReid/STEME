/** Copyright John Reid 2011, 2012
 *
 * \file Exposes background parts of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/background_model.h>
#include <steme/background_model_on_the_fly.h>
#include <steme/data.h>

#include <myrrh/python/boost_range.h>

#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace py = boost::python;
using namespace steme;
using namespace steme::python;

typedef data< STEME_INDEX_MODULE_SPEC > data_t;
typedef background_model_meta< STEME_INDEX_MODULE_SPEC >::type bg_model_t;


namespace steme {
const double log_quarter = std::log( 0.25 ); /**< log(.25) */
} //namespace steme



double
get_freq( const zero_order_frequencies & freqs, size_t base ) {
	if( base >= freqs.dist.size() ) {
		throw std::out_of_range( "Index into frequencies out of range." );
	}
	return freqs.dist[ base ];
}


zero_order_frequencies::ptr
zero_order_from_seq( py::object py_seq )
{
    if( 4 != boost::python::len( py_seq ) ) {
        throw std::logic_error( "Wrong number of elements in sequence." );
    }
    return zero_order_frequencies::ptr(
        new zero_order_frequencies(
            myrrh::python::make_extract_iterator< double >( myrrh::python::py_seq_iterator( py_seq, 0 ) ),
            myrrh::python::make_extract_iterator< double >( myrrh::python::py_seq_iterator( py_seq, 4 ) )
        )
    );
}




#if STEME_INDEX_MODULE_TYPE == STEME_INDEX_DEFAULT_TYPE


bg_model::ptr
create_bg_model_from_0_order(
    size_t W,
    data_t & data,
    const zero_order_frequencies & freqs
)
{
    bg_model::ptr result( new bg_model( W, freqs) );
    markov_0_order_likelihoods_calculator( data, freqs, result->wmer_LLs );
    return result;
}




bg_model::ptr
create_bg_model_from_base_likelihoods(
    size_t W,
    data_t & data,
    const likelihoods_vec_vec_t * base_likelihoods,
    const zero_order_frequencies & freqs
)
{
    bg_model::ptr result( new bg_model( W, freqs, base_likelihoods ) );
    base_to_wmer_likelihoods_calculator( data, *result->base_LLs, result->wmer_LLs );
    return result;
}



template< size_t Order >
bg_model::ptr
create_bg_model_from_Markov_model(
    size_t W,
    data_t & data,
    const complete_markov_model< Order, 4, double > & model,
    const zero_order_frequencies & freqs
)
{
    bg_model::ptr result( new bg_model( W, freqs ) );
    markov_wmer_likelihoods_calculator< complete_markov_model< Order, 4, double > >( data, model, result->wmer_LLs );
    return result;
}



#else //STEME_INDEX_MODULE_TYPE == STEME_INDEX_DEFAULT_TYPE

template< typename MultiArray >
double
get_smallest_ll( const MultiArray & x ) {
    if( x.empty() ) {
        throw std::logic_error( "Markov model has no likelihoods." );
    }
    return *std::min_element( x.data(), x.data() + x.num_elements() );
}



/// Calculates the likelihood of a infix in context on the fly.
template< typename MarkovModel >
struct markov_on_the_fly_evaluator
: on_the_fly_evaluator< STEME_INDEX_MODULE_SPEC >
{
    const MarkovModel & mm; ///< The Markov model
    double smallest_ll; ///< The smallest log likelihood in the Markov model.

    markov_on_the_fly_evaluator( const MarkovModel & mm )
    : mm( mm )
    , smallest_ll( get_smallest_ll( mm.mm.storage ) )
    { }

    virtual
    ~markov_on_the_fly_evaluator()
    { }

    struct addition
    {
        addition( double & x ) : x( x ) { }
        void operator()( double y ) const { x += y; }
        double & x;
    };

    virtual
    double
    evaluate( it context, it begin, it end ) {
        double result = 0.;
        mm.evaluate_in_context(
            context,
            begin,
            end,
            boost::make_function_output_iterator( addition( result ) )
        );
        return result;

    }

    virtual
    double
    min_ll( size_t w, size_t W ) {
        return ( W - w ) * smallest_ll;
    }
};




template< size_t Order >
on_the_fly_bg_model< STEME_INDEX_MODULE_SPEC >::ptr
create_on_the_fly_bg_model_from_Markov_model(
    const complete_markov_model< Order, 4, double > & model
) {
    typename on_the_fly_evaluator< STEME_INDEX_MODULE_SPEC >::ptr evaluator( new markov_on_the_fly_evaluator< complete_markov_model< Order, 4, double > >( model ) );
    return on_the_fly_bg_model< STEME_INDEX_MODULE_SPEC >::ptr( new on_the_fly_bg_model< STEME_INDEX_MODULE_SPEC >( evaluator ) );
}


#endif //STEME_INDEX_MODULE_TYPE == STEME_INDEX_DEFAULT_TYPE




void
expose_bg() {

#ifdef STEME_USE_OLD_BG_MODEL
	py::scope().attr("_using_old_bg_model") = true;
	std::cout << "WARNING: Using old background model. If you did not intend this then check your configuration!" << std::endl;
#else //STEME_USE_OLD_BG_MODEL
	py::scope().attr("_using_old_bg_model") = false;
#endif //STEME_USE_OLD_BG_MODEL


	//
	// zero_order_frequencies
	//
	py::class_< zero_order_frequencies, zero_order_frequencies::ptr > zero_order_freqs_class(
		"ZeroOrderFrequencies",
		"0-order frequencies.",
		py::no_init
    );
	zero_order_freqs_class.def( "__init__", make_constructor( zero_order_from_seq ), "Constructs from the 0-order occurrences in the sequence." );
	zero_order_freqs_class.def( "add_pseudo_counts", &zero_order_frequencies::add_pseudo_counts, "Add pseudo-counts and return a new 0-order frequencies object." );
	zero_order_freqs_class.def( "freq", get_freq, "Get the frequency of the base." );
	zero_order_freqs_class.def_readonly( "total_counts", &zero_order_frequencies::total_counts, "Total counts used to make this frequencies object." );





    //
    // Background model
    //
    py::class_<
        bg_model_t,
        boost::noncopyable,
        bg_model_t::ptr
    > bg_model_class(
        "BgModel",
        "Models the likelihoods of W-mers under a background model.",
        py::no_init
    );
#if STEME_INDEX_MODULE_TYPE == STEME_INDEX_DEFAULT_TYPE
    bg_model_class.def_readonly( "freqs", &bg_model_t::freqs, "The background frequencies." );


    py::def(
        "create_bg_model_from_0_order",
        create_bg_model_from_0_order,
        (
            py::arg( "W" ),
            py::arg( "data" ),
            py::arg( "freqs" )
        ),
        "Creates a background model from 0-order frequencies."
    );



    py::def(
        "create_bg_model_from_base_likelihoods",
        create_bg_model_from_base_likelihoods,
        (
            py::arg( "W" ),
            py::arg( "data" ),
            py::arg( "base_LLs" ),
            py::arg( "freqs" )
        ),
        "Creates a background model from likelihoods given per-base.",
        py::with_custodian_and_ward_postcall<
            0,
            2,
            py::with_custodian_and_ward_postcall< 0, 3 >
        >() // don't let base_LLs or data be destroyed before the background model
    );


/// Registers a create background from Markov model of given order
# define REGISTER_CREATE_BG_FROM_MARKOV_MODEL_OF_ORDER( z, Order, text ) \
    py::def( \
        MYRRH_MAKE_STRING( "create_bg_model_from_Markov_model_" << Order ).c_str(), \
        create_bg_model_from_Markov_model< Order >, \
        ( \
            py::arg( "W" ), \
            py::arg( "data" ), \
            py::arg( "model" ), \
            py::arg( "freqs" ) \
        ), \
        "Creates a background model from the Markov model." \
    );

    // define create from Markov model functions.
    BOOST_PP_REPEAT_FROM_TO( 0, STEME_MAX_MARKOV_MODEL_ORDER, REGISTER_CREATE_BG_FROM_MARKOV_MODEL_OF_ORDER, )




#else //STEME_INDEX_MODULE_TYPE == STEME_INDEX_DEFAULT_TYPE


/// Registers a create background from Markov model of given order
# define REGISTER_CREATE_BG_FROM_MARKOV_MODEL_OF_ORDER( z, Order, text ) \
    py::def( \
        MYRRH_MAKE_STRING( "create_bg_model_from_Markov_model_" << Order ).c_str(), \
        create_on_the_fly_bg_model_from_Markov_model< Order >, \
        ( \
            py::arg( "model" ) \
        ), \
        "Creates a background model from the Markov model." \
    );


    // define create from Markov model functions.
    BOOST_PP_REPEAT_FROM_TO( 0, STEME_MAX_MARKOV_MODEL_ORDER, REGISTER_CREATE_BG_FROM_MARKOV_MODEL_OF_ORDER, )


#endif //STEME_INDEX_MODULE_TYPE == STEME_INDEX_DEFAULT_TYPE

}
