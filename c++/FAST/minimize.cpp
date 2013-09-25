#include "minimize.h"
#include "mnbrak.h"
#include "brent.h"

//
// Value is set to argmin_x f(x) and f(value) is returned.
//
double
minimize(
    double
    (*f)(
        double ), double& value ) {

    //
    // Function minimization using "Numerical Recipes" functions.
    //

    double lower_bnd = -1;
    double middle = 0.5;
    double upper_bnd = 3.0;
    double tolerance = 1e-6;
    double fa, fb, fc;

    mnbrak(
        &lower_bnd, &middle, &upper_bnd, &fa, &fb, &fc, f );
    return brent(
        lower_bnd, middle, upper_bnd, f, tolerance, &value );
}
