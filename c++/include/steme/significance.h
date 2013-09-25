/** Copyright John Reid 2011
 *
 * \file
 * \brief Defines the model type for STEME algorithm.
 *
 */

/**
 * Some of the code in this file was derived from MEME's source.
 * Here is the copyright notice:

    Copyright  (c)  1994-2009  The  Regents  of the  University of
    California.  All  Rights  Reserved.

    Permission  to use,  copy,  modify,  and  distribute  any part
    of this  software for  educational,  research  and  non-profit
    purposes,  without  fee,  and  without a written  agreement is
    hereby  granted,  provided  that the  above  copyright notice,
    this paragraph  and the following  three  paragraphs appear in
    all copies.

    Those  desiring to  incorporate this  software into commercial
    products  or use for  commercial  purposes  should contact the
    Technology  Transfer  Office,  University of California,   San
    Diego,  9500 Gilman Drive,  La Jolla,  California, 92093-0910,
    Phone: (858) 534-5815.

    IN  NO  EVENT  SHALL THE  UNIVERSITY  OF CALIFORNIA  BE LIABLE
    TO  ANY  PARTY FOR  DIRECT,  INDIRECT, SPECIAL, INCIDENTAL, OR
    CONSEQUENTIAL  DAMAGES,  INCLUDING  LOST PROFITS, ARISING  OUT
    OF  THE  USE  OF  THIS  SOFTWARE,  EVEN  IF THE UNIVERSITY  OF
    CALIFORNIA  HAS  BEEN  ADVISED  OF  THE  POSSIBILITY  OF  SUCH
    DAMAGE.

    THE SOFTWARE  PROVIDED HEREUNDER IS ON AN  "AS IS" BASIS,  AND
    THE  UNIVERSITY OF CALIFORNIA  HAS  NO OBLIGATIONS  TO PROVIDE
    MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
    THE UNIVERSITY  OF CALIFORNIA  MAKES  NO  REPRESENTATIONS  AND
    EXTENDS  NO  WARRANTIES  OF  ANY  KIND,  EITHER  EXPRESSED  OR
    IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
    OF  MERCHANTABILITY  OR  FITNESS FOR A  PARTICULAR PURPOSE, OR
    THAT  THE USE  OF THE MATERIAL  WILL NOT  INFRINGE ANY PATENT,
    TRADEMARK OR OTHER RIGHTS.

 */


#ifndef STEME_JR_31AUG2011_SIGNIFICANCE_H_
#define STEME_JR_31AUG2011_SIGNIFICANCE_H_

#include <steme/model.h>

#include <pvalues/calculator.h>

namespace steme {



/**
 * Measures a model's significance.
 */
template< typename Spec = default_spec >
struct significance
: boost::noncopyable
{
    typedef significance< Spec >                   self_t;                 ///< This type.
    typedef boost::shared_ptr< self_t >            ptr;                    ///< Shared pointer to this type.
    typedef data< Spec >                           data_t;                 ///< The data type.
    typedef model< Spec >                          model_t;                ///< The model type.

    data_t &                                       _data;                  ///< The data.
    zero_order_frequencies &                       _freqs;                 ///< The 0-order frequencies.
    pvalues::llr_pvalue_calculator::shared_ptr     _calculator;            ///< Calculates p-values.

    /// Constructor.
    significance(
        data_t & _data,
        zero_order_frequencies & _freqs,
        pvalues::llr_pvalue_calculator::shared_ptr _calculator
    )
    : _data( _data )
    , _freqs( _freqs )
    , _calculator( _calculator )
    { }


    /**
     * The logarithm of the E-value of the model.
     */
    double
    log_E_value( const model_t & model ) {
        const double log_pop = log_product_p_values( model );
        // pretend that we only had the maximum we can calculate significance for if we actually have too many
        const size_t num_samples = std::min( _calculator->get_max_N(), size_t( model.bs.pssm.num_samples ) );
        return pvalues::log_qfast( model.W(), log_pop ) + pvalues::get_log_nalign( num_samples, _data.get_occurrence_count( model.W() ) );
    }


    /**
     * The log product of p-values of the model (as in MEME).
     */
    double
    log_product_p_values( const model_t & model ) {
        MYRRH_ASSERT( log_probs_sum_to_1( _freqs.dist_logs ) );
        double log_pop = 0.;
        // pretend that we only had the maximum we can calculate significance for if we actually have too many
        const size_t num_samples = std::min( _calculator->get_max_N(), size_t( model.bs.pssm.num_samples ) );
        BOOST_FOREACH( PssmStorage::matrix::const_reference pssm_col, model.bs.pssm.log_probs.m ) {
            const double ic = calculate_column_ic( pssm_col, _freqs.dist_logs.data() );
            const double llr = ic * num_samples;
            const double log_pv = ( *_calculator )( num_samples, llr );
            log_pop += log_pv;
        }
        return log_pop;
    }


};


} // namespace steme

#endif // STEME_JR_31AUG2011_SIGNIFICANCE_H_
