#ifndef _ENTROPY_H_
#define _ENTROPY_H_

#define NOT_COMPUTED 0

//
// Function: Computes the lattice pmf (log-values in an array of size
//           L*Q) for a sum of scores based on the sFFT algorithm.
//
// Parameters: log_pmf = Log of lattice pmf for single sample
//             L = Number of i.i.d. samples
//             s = Score for which theta is optimized
//             Q = Lattice size for single sample
//             step = Step size for the lattice
//             theta = Optimal theta shift
//             log_mgf = Log of the mgf of the pmf at theta
//
// If theta is set to NOT_COMPUTED then it is computed and otherwise
// the given theta is used.
//
double*
sFFT(
    double* log_pmf, int L, double s, int Q, double step, double& theta,
    double& log_mgf );

//
// Function: Computes the lattice pmf (log-values in an array of size
//           L*Q) for a sum of scores based on the csFFT algorithm.
//
// Parameters: log_pmf = Log of lattice pmf for single sample
//             L = Number of i.i.d. samples
//             s = Score for which theta is optimized
//             Q = Lattice size for single sample
//             step = Step size for the lattice
//             theta = Optimal theta shift
//             log_mgf = Log of the mgf of the pmf at theta
//
// If theta is set to NOT_COMPUTED then it is computed and otherwise
// the given theta is used.
//
double*
csFFT(
    double* log_pmf, int L, double s, int Q, double step, double& theta,
    double& log_mgf );

//
// Function: Computes the p-value for a sum of scores based on the lattice-refinement 
//           version of the csFFT algorithm.
//
// Parameters: log_pmf = Log of lattice pmf for single sample
//             L = Number of i.i.d. samples
//             s = Score for which theta is optimized
//             Q = Lattice size for single sample
//             step = Step size for the lattice
//             accuracy = Desired number of decimal places of accuracy for the p-value
//             theta = Optimal theta shift
//             log_mgf = Log of the mgf of the pmf at theta
//             log_max_pval = Log of an upper bound on the p-value
//             log_min_pval = Log of a lower bound on the p-value
//
// If theta is set to NOT_COMPUTED then it is computed and otherwise
// the given theta is used.
//
void
refine_csFFT(
    double* log_pmf, int L, double s, int Q, double step, double accuracy,
    double& theta, double& log_mgf, double& log_max_pval, double& log_min_pval );

#endif
