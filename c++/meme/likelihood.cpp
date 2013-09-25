/**
 * The code in this file was derived from MEME's source.
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


#include <cmath>
#include <limits>

/**********************************************************************/
/*
	log_qfast

	Calculate the log p-value of the log of the
	product of uniform [0,1] random variables.

*/
/**********************************************************************/
double
log_qfast(
	int n,               /* number of random variables in product */
	double logk          /* product of random variables */
) {
	int i;
	double term, phi;

	if (n == 0)
		return 0; /* worst possible log p-value */

	phi = term = 1;
	for( i = 1; i != n; i++ ) {
		term *= -logk / i;
		phi += term;
	}

	return logk + std::log( phi );
}                               /* qfast */



/**********************************************************************/
/*
	get_log_nalign

	Get an upper bound on the number of independent alignments
	of segments of length w.

	Same as get_log_nalign() in meme source when model is tcm and invcomp is true.
*/
/**********************************************************************/
double
get_log_nalign(
	int N,                      /* number of occurrences */
	int num_possible_sites      /* number of possible sites */
)
{
	double log_nalign = 0.;        /* log number alignments */

    if( N > num_possible_sites ) {                       /* impossible N */

    	log_nalign = std::numeric_limits< double >::max();

    } else { /* remove 1 site per site */

    	for( int i = 0; i != N; ++i ) {
			log_nalign += std::log( (num_possible_sites - i) * 2 / (i + 1) );
		}

    }

	return log_nalign;
}                               /* double get_log_nalign */



extern "C" {

/**********************************************************************/
/*
	get_log_sig

  	Calculate the statistical significance of the alignment given
	its score and the type of objective function in force.

	If N>0, returns log E-value.
	If N==0, returns log p-value.
*/
/**********************************************************************/
double
get_log_sig (
	double score,          /* score of alignment */
	int w,                 /* width of motif */
	int N,                 /* number of sites */
	int num_possible_sites /* number of possible sites */
  )
{
	double log_pv = log_qfast( w, -score );       /* p-value of product of p-values */

	if( N ) {                           /* use E-value of alignment */
		return log_pv + get_log_nalign( N, num_possible_sites );
	} else {                           /* use p-value of alignment */
		return log_pv;
	}
}                               /* get_log_sig */

}
