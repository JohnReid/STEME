/**
 * Copyright John Reid 2011, 2012
 *
 * \file
 * \brief Code to test implementation of Bejerano's convex algorithm for calculating p-values.
 *
 */

#define BOOST_TEST_MODULE Bejerano convex pvalues test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/timer.hpp>

#include <pvalues/bejerano.h>

#include <iostream>

using namespace pvalues;
using namespace std;


BOOST_AUTO_TEST_CASE( binary_search_1 )
{
    using namespace boost::lambda;
    const int min = 20;
    for( int max = min; min + 5 != max; ++max ) {
        std::cout << "[" << min << "," << max << ")\n";
        for( int threshold = min - 1; max + 2 != threshold; ++threshold ) {
            const int result = do_binary_search( _1 >= threshold, min, max );
            std::cout << threshold << " : " << result << "\n";
            BOOST_CHECK(
                ( threshold < min && result == min )
                || ( max < threshold && result == max )
                || result == threshold
            );
        }
    }
}




namespace pvalues {
namespace spec {

/// Specification for our testing.
struct Test { };

/** Float meta-function specialisation. */
template<>
struct Float< Test > {
    typedef double type;
};

/** Log on/off meta-function specialisation. */
template<>
struct LogOn< Test > {
    static const bool value = false;
};

} // namespace pvalues
} // namespace spec

typedef bejerano< pvalues::spec::Test > pval_calculator_t;
//typedef bejerano<> pval_calculator_t;



template< typename QRange >
void
check_p_value( unsigned n, double llr, const QRange & q, double target_pvalue, double percent_tol = 1e-4 ) {

    cout << "\n";
    cout << "n   = " << n   << "\n";
    cout << "llr = " << llr << "\n";

    pval_calculator_t::log_arithmetic_t log_arith( n + 1 );
    pval_calculator_t pval_calculator( log_arith, n, llr, q );

    boost::timer t;
    const pval_calculator_t::float_t log_pval = pval_calculator();
    cout << "Took " << t.elapsed() << " seconds.\n";
    const double pval = log_arith.exp( log_pval );
    const double relative_error = std::fabs( pval - target_pvalue ) / std::max( target_pvalue, pval );

    cout << "target p-value = " << target_pvalue  << "\n";
    cout << "log(p-value)   = " << log_pval       << "\n";
    cout << "p-value        = " << pval           << "\n";
    cout << "relative error = " << relative_error << "\n";
    cout << "\n";

    // check is within 1e-4% of target value
    BOOST_CHECK_CLOSE( pval, target_pvalue, percent_tol );
}



// test we can handle cases where LLR = 0
BOOST_AUTO_TEST_CASE( bejerano_0_LLR )
{
    const double q[] = { .2, .2, .3, .3 };
    const unsigned n = 10;
    const double target_pvalue = 0.970607; // shouldn't this be 1? Bejerano's implementation gives this value.

    check_p_value( n, 0., q, target_pvalue );
}


/// \todo Add example where some of the intermediary values are integers. This would test
/// corner cases for rounding in algorithm

// example from paper:
//   Efficient Exact p-value Computation
//   and Applications to Biosequence Analysis
BOOST_AUTO_TEST_CASE( bejerano_paper_example )
{
    const double q[] = { .1, .45, .45 };
    const unsigned n = 2;
    const double llr = 3.4 / 2;
    const double target_pvalue = .19;

    check_p_value( n, llr, q, target_pvalue );
}


BOOST_AUTO_TEST_CASE( bejerano_simple )
{
    const double q[] = { .2, .2, .3, .3 };
    const unsigned n = 10;
    const double llr = 1.4 * n;
    const double target_pvalue = 2.048e-07;

    check_p_value( n, llr, q, target_pvalue );
}

#ifdef NDEBUG
BOOST_AUTO_TEST_CASE( bejerano_longer )
{
    std::cout << "sizeof(float) = " << sizeof(float) << "\n";
    std::cout << "sizeof(double) = " << sizeof(double) << "\n";
    std::cout << "sizeof(long double) = " << sizeof(long double) << "\n";

    const double q[] = { .2, .2, .3, .3 };
    const unsigned n = 1200;
    const double llr = 600.;
    const double target_pvalue = 8.07923e-260;

    check_p_value( n, llr, q, target_pvalue, 1. ); // use 1% tolerance
}
#endif //NDEBUG


