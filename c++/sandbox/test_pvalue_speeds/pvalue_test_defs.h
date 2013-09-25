/**
 * Copyright John Reid 2011
 *
 * \file Common definitions for p-value testing
 */



#ifndef SANDBOX_TEST_PVALUE_SPEEDS_DEFS_H
#define SANDBOX_TEST_PVALUE_SPEEDS_DEFS_H

#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

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
read_p_value_test_cases( std::ifstream & in, p_value_test_case_vec & test_cases, size_t max_N = 0, size_t max_cases = 0 );



/**
 * A function that calculates a p-value for a test case.
 */
struct pvalue_calculator : std::unary_function< const p_value_test_case &, double >
{
	typedef boost::shared_ptr< pvalue_calculator > ptr;
	virtual double operator()( const p_value_test_case & args ) = 0;
};



pvalue_calculator::ptr
create_PVAL_pvalue_calculator(
	double * pu,
	int alg,
	bool use_qfast
);

pvalue_calculator::ptr
create_MEME_pvalue_calculator(
	int minN,
	int maxN,
	double * pu,
	int range,
	bool use_qfast
);

pvalue_calculator::ptr
create_FAST_pvalue_calculator(
	int minL,
	int maxL,
	double * pu,
	int Q = 16384
);



void
init_MEME_llr_pv_tables(
  int min,				              /* minimum number of sites */
  int max,				              /* maximum number of sites */
  double * pu                         /* Background frequencies. */
);









#endif //SANDBOX_TEST_PVALUE_SPEEDS_DEFS_H

