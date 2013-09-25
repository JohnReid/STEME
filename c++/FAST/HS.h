#ifndef _HS_H_
#define _HS_H_

//
// Function: Computes and returns the entire lattice pmf (log-values in an array of
//           size Q) for the multinomial llr score based on the shifted-hirji algorithm
//           using lists.
//
// Parameters: N = Number of objects 
//             K = Number of bins 
//             pu = Null multinomial distribution 
//             Q = Lattice size 
//             step = Step size for the lattice (set by the function)
//             index = Array of indices for non-zero pmf entries (set by the function)
//             size = Size of the index array (set by the function)
//
double*
HSlist(
    int N, int K, double *pu, int Q, double& step, int*& index, int& size );

#endif 
