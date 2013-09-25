/**
 * Copyright John Reid 2011, 2013
 *
 * \file Common definitions for p-value testing
 */



#ifndef STEME_PVALUE_JR_15AUG2011_DEFS_H
#define STEME_PVALUE_JR_15AUG2011_DEFS_H

#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>

#include <cmath>
#include <functional>
#include <numeric>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>


//
// A test case for p-value generation
//
struct p_value_test_case {
	typedef std::vector< double > IC_vec;

	size_t N;
	IC_vec ICs;

	size_t W() const { return ICs.size(); }
	double IC() const { return std::accumulate( ICs.begin(), ICs.end(), 0.0 ); }

	p_value_test_case( size_t N=0 ) : N( N ) { }
};

/// Stream test cases
std::ostream &
operator<<( std::ostream & os, const p_value_test_case & args );

/// Vector of test cases
typedef std::vector< p_value_test_case > p_value_test_case_vec;

/// Read test cases from file.
void
read_p_value_test_cases(
    std::ifstream & in,
    p_value_test_case_vec & test_cases,
    size_t min_N = 0,
    size_t max_N = 0,
    size_t max_cases = 0
);



/**
 * A function that calculates a p-value for a test case.
 */
struct pvalue_calculator : std::unary_function< const p_value_test_case &, double >
{
	virtual ~pvalue_calculator() {}

	typedef boost::shared_ptr< pvalue_calculator > ptr;
	virtual double operator()( const p_value_test_case & args ) = 0;
};



pvalue_calculator::ptr
create_PVAL_pvalue_calculator(
	const double * pu,
	int alg,
	bool use_qfast
);

pvalue_calculator::ptr
create_MEME_pvalue_calculator(
	int minN,
	int maxN,
	const double * pu,
	int range,
	bool use_qfast
);

pvalue_calculator::ptr
create_HS_pvalue_calculator(
    int maxN,
    const double * pu,
    size_t Q,
    bool use_qfast
);

pvalue_calculator::ptr
create_FAST_pvalue_calculator(
	const double * pu,
	int Q = 16384,
    bool use_qfast = false
);



void
init_MEME_llr_pv_tables(
  int min,				              /* minimum number of sites */
  int max,				              /* maximum number of sites */
  double * pu                         /* Background frequencies. */
);











#endif //STEME_PVALUE_JR_15AUG2011_DEFS_H

