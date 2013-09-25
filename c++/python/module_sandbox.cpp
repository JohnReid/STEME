/** Copyright John Reid 2009, 2012
 *
 * \file Sandbox python extension.
 */

#include <boost/python.hpp>

namespace py = boost::python;

struct held {
	typedef boost::shared_ptr< held > ptr;

	std::string _data;
};

struct holder {
	held _held;

	held & get_held() { return _held; }
};




BOOST_PYTHON_MODULE( _sandbox )
{
	py::class_< held, held::ptr > held_class(
		"Held",
		"Class held in another.",
		py::no_init
	);
	held_class.def_readwrite( "data", &held::_data, "Data." );

	py::class_< holder > holder_class(
		"Holder",
		"Class that holds some data."
	);
	//holder_class.def_readonly( "held", &holder::get_held, "Held object." );
	holder_class.def( "held_by_value", &holder::get_held, "Held object.", py::return_value_policy< py::return_by_value >() );
	holder_class.def( "held_internal_ref", &holder::get_held, "Held object.", py::return_internal_reference<>() );
	holder_class.add_property( "held", py::make_getter( &holder::_held ), "Held object." );
}
