/**
 * Copyright John Reid 2011
 */

#include "pvalue_test_defs.h"

using namespace std;

ostream &
operator<<( ostream & os, const p_value_test_case & args ) {
    os << "N=" << args.N << "; W=" << args.W() << "; IC=";
    std::copy( args.ICs.begin(), args.ICs.end(), std::ostream_iterator< double >( os, "," ) );
    return os;
}

