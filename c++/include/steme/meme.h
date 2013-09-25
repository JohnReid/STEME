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


#ifndef STEME_JR_13AUG2011_MEME_H_
#define STEME_JR_13AUG2011_MEME_H_

#define MSN	 24 		/* maximum length of sample name */
				/* MSN + 40 < PAGEWIDTH (see meme.h) */
#define MAXALPH  28		/* maximum length of alphabet + 1 for 'X' */
#define MAXG 101		/* maximum number of motifs + 1 */
#define MAXSITE 300		/* maximum length of a site */
#define MINSITES 2		/* minimum number of sites in valid motif */
#define LLR_RANGE 200		/* range of scaled LLR statistic */

#define MINCONS 0.2		/* Display 'X' as consensus if no letter f > */
#define LOGOHEIGHT 7.5		// height of sequence logo in cm.
#define MAXLOGOWIDTH 30		// maximum width of sequence logo in cm.

/* minimum allowable motif width before shortening;
   never make less than 2 or will crash! */
#define MIN_W 8
/* maximum allowable length before shortening */
#define MAX_W 50

#define MNAME 20		/* names of known motifs */
#define NMOTIFS MAXG		/* maximum number of known motifs */


/* default size of heap for branching search */
#define HSIZE 64
#define HS_DECREASE 2

/* default branching factor for branching search */
#define BFACTOR 3

/* Amount of error tolerated in probability column sums (they should sum to
   approximately 1): */
#define ERR_EPSILON 0.01

/*
  round x to d significant digits and put in y
*/
#define RNDDIG 14
#define RND(x, d, y) {                                                  \
  if (x > 0) {                                                          \
    double _z_ = exp10(ceil((d)-1-log10(x)));                           \
    y = rint(_z_*(x))/_z_;                                              \
  } else if (x < 0) {                                                   \
    double _z_ = exp10(ceil((d)-1-log10(-x)));                          \
    y = -rint(_z_*(-x))/_z_;                                            \
  } else {                                                              \
    y = 0;                                                              \
  }                                                                     \
}

typedef int BOOLEAN;

/* function prototypes */
extern void init_log (void);
extern void init_exp (void);

extern void reset_llr_pv_tables();

extern
void
init_llr_pv_tables(
  int min,				/* minimum number of sites */
  int max,				/* maximum number of sites */
  int alength,				/* alphabet length */
  double *back,				/* background frequencies */
  BOOLEAN pal				/* sites are palindromes */
);

extern
double
get_llr_pv(
  double llr,				/* log likelihood ratio */
  double n,				/* number of sequences in alignment */
  int w,				/* width of alignment */
  int range,				/* desired range of scaled LLR */
  double frac,				/* speedup factor */
  int alength,				/* length of alphabet */
  double *dd 				/* alphabet frequency distribution */
);

extern
double
get_llr_mean(
  double n 				/* number sequences in alignment */
);

/// Get the log significance. MEME function.
extern
double
get_log_sig (
	double score,          /* score of alignment */
	int w,                 /* width of motif */
	int N,                 /* number of sites */
	int num_possible_sites /* number of possible sites */
);


#endif // STEME_JR_13AUG2011_MEME_H_

