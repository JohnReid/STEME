#ifndef _MINIMIZE_H_
#define _MINIMIZE_H_

//
// Value is set to argmin_x f(x) and f(value) is returned.
//
double
minimize(
    double
    (*f)(
        double ), double& value );

#endif
