/** Copyright John Reid 2011
 *
 * \file Exposes MEME parts of STEME algorithm to python.
 */

#include "steme_python_defs.h"

#include <steme/background_model.h>

extern "C" {
#include "steme/meme.h"
}

namespace py = boost::python;
using namespace steme;
using namespace steme::python;


namespace steme {

/// The largest number of samples we are prepared to do an E-value calculation for.
size_t max_samples_in_E_value_calc = 4000;

/** Stores min and max values of llr p-value tables. */
boost::optional< boost::tuple< int, int > >
MEME_llr_pv_tables_initialised;

} // namespace steme




/**
 * Reset MEME LLR pv tables. Will cause memory leak in MEME code.
 */
void
reset_MEME_llr_pv_tables() {
	MEME_llr_pv_tables_initialised = boost::optional< boost::tuple< int, int > >();
	reset_llr_pv_tables();
}



/**
 * Initialise MEME LLR pv tables.
 */
void
init_MEME_llr_pv_tables(
	int min,                              /* minimum number of sites */
	int max,                              /* maximum number of sites */
	zero_order_frequencies & bg           /* Background frequencies. */
) {
	//check not initialised already.
	if( MEME_llr_pv_tables_initialised ) {
		throw std::logic_error( "MEME LLR p-value tables already initialised." );
	}

	//check we aren't trying to do something that will take too long.
	if( max > int( max_samples_in_E_value_calc ) ) {
		throw std::logic_error(
			MS_MAKE_STRING(
				"Cannot calculate E-value for models with number of samples ("<<max<<") >"<<max_samples_in_E_value_calc
			)
		);
	}

	//initialise the tables.
	MEME_llr_pv_tables_initialised = boost::make_tuple( min, max );
	init_log();
	init_exp();
	init_llr_pv_tables( min, max, 4, bg.dist.data(), 0 );
}



void
expose_meme() {

	py::def( "reset_MEME_llr_pv_tables", reset_MEME_llr_pv_tables, "Reset MEME LLR p-value tables." );
	py::def( "init_MEME_llr_pv_tables", init_MEME_llr_pv_tables, "Initialise MEME LLR p-value tables." );

}
