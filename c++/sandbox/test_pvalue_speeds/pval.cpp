/*-----------------------------------------------------------------------------

 G^2 p-value computation

 Accompanies application note submission:
 Gill Bejerano, "Branch and Bound Computation of Exact P-Values",
 Bioinformatics, submitted, 2006.

 Implements methods described in detail in:
 G. Bejerano, N. Friedman and N. Tishby. "Efficient exact p-value computation
 for small sample, sparse and surprising categorical data", J. Computational
 Biology, 11(5): 867-886, 2004.

 Please quote the JCB paper if you find this code/approach useful.

 See http://www.soe.ucsc.edu/~jill/ for additional material.

 -----------------------------------------------------------------------------*/

#include "pvalue_test_defs.h"

#include <boost/foreach.hpp>

#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

using namespace std;

#define PRINT_DOUBLE false
// toggles stdout printout - double for accuracy or %f for readability

#define FORBID_SUBTRACTIONS true
// allowing subtractions is equivalent theoretically,
// speeds up computation, but sometimes introduces numerical errors

#define F(X) X*(lg[X]-lgQ[i])+(n_comp-X)*(lg[n_comp-X]-lgq_star[i])

// globals save time passing them between recursion calls

#ifdef LOOKUP
// you can add this flag at compilation time (see README)
const int D_LOGADD = 30; // y/x < 1e-D discard
const int G_LOGADD = 100; // step function look-up every 1/G
const int D_LOGSUB = 17; // above 16.25.. the func is zero
const int G_LOGSUB = 200;

int d_logadd = int(D_LOGADD*log(10)*G_LOGADD);
int d_logsub = int(D_LOGSUB*log(10)*G_LOGSUB);
double *logaddf,*logsubf;
#endif



double d_min(double *d, int from, int to) {
	double min = d[from];
	for (int i = from + 1; i <= to; i++)
		if (d[i] < min)
			min = d[i];
	return (min);
}

double chi2pval(int df, double value) {
	double gammq(double a, double x);

	if (value < 0.0) {
		cerr << "Warning: in chi2pval rounding " << value << " to zero."
				<< endl;
		value = 0.0;
	}
	return (gammq(df / 2.0, value / 2.0));
}

// inline double nchoosek (int n, int k)
// {
//   double result = 1.0;
//   if (n-k < k)
//     k = n-k;
//   for (int i = 0; i < k; i++)
//     result *= double(n-i)/(k-i);
//   return(result);
// }


// long nodes_rec(int n, int k)
// //very inefficient, opens up the whole tree...
// {
//   if (k<2)
//     return (k);
//   long sum = n+1;
//   for (int i=0; i<=n; i++)
//     sum += nodes_rec(i, k-1);
//   return(sum);
// }


char * currtime(void) {
	time_t now = time(NULL);
	return (ctime(&now));
}



//int compare_forperm(const void *i, const void *j) {
//	double diff = Q0[*(int *) i] - Q0[*(int *) j];
//	return (diff < 0.0 ? -1 : (diff > 0.0) ? +1 : 0);
//}

/**
 * A wrapper around Bejerano's pvalue calculations
 */
struct pval_pvalue_calculator : pvalue_calculator {
	int alg;

	int N, K, *C0, *n, *perm, underflow, biggest_N;
	double *Q0, *lgQ0, eps, KL_thresh, lg_pval, pval;
	double *lg, *cumlg, *Q, *lgQ, *lgQcomp, *lgQcompmin, simsec, runtime;
	double *Qcomp, *Qcompmin;
	long reccalls, logops, discards;
	bool *permver;
	bool use_qfast;

	pval_pvalue_calculator( int alg, double * pu, bool use_qfast )
	: alg(alg)
	, use_qfast( use_qfast )
	{
		switch( alg ) {
		case 0: cout << "Using chi^2 approximation.\n"; break;
		case 1: cout << "Using exact method.\n"; break;
		case 2: cout << "Using branch-and-bound method.\n"; break;
		case 3: cout << "Using convex method.\n"; break;
		default: cout << "Unknown method!\n"; break;
		}

		if( use_qfast ) cout << "Using QFAST to combine column p-values.\n";
		else cout << "Not using QFAST to combine column p-values.\n";

		int input, counter = 0, df;
		unsigned int seed;

		cerr.precision(10);

		seed = (unsigned int) time(NULL);
		cout << "seed = " << seed << "\n";
		srand(seed);

		K = 4;

		Q0 = new double[K + 1]; // H_0: Q0[1], .., Q0[K]
		lgQ0 = new double[K + 1]; // mainly for simulation

		C0 = new int[K + 1]; // observed type, w.r.t Q0
		n = new int[K + 1]; // recursive partial type for 'expand'
		// coord[0] unused

		perm = new int[K + 1];
		Q = new double[K + 1]; // permuted H_0
		lgQ = new double[K + 1];
		lgQcomp = new double[K + 1]; // lgQcomp[i] = log(Q[i+1]+...+Q[K])
		lgQcompmin = new double[K + 1]; // lgQcompmin[i] = log(min(Q[i+1],..,Q[K]))
		Qcomp = new double[K + 1]; // for convex binsearch
		Qcompmin = new double[K + 1]; // for convex binsearch
		// entries coord[0], lgQcomp[K], lgQcompmin[K] are not used

	#ifdef LOOKUP
		cout << "Creating log lookup tables.\n";

		logaddf = new double [d_logadd+1];
		for (int i=0; i<= d_logadd; i++)
		logaddf[i] = log(1+exp(-double(i)/G_LOGADD));

		logsubf = new double [d_logsub+1];
		logsubf[0] = -HUGE_VAL;
		for (int i=1; i<= d_logsub; i++)
		logsubf[i] = log(1-exp(-double(i)/G_LOGSUB));
	#else //LOOKUP
		cout << "Not using log lookup tables.\n";
	#endif

		if (FORBID_SUBTRACTIONS)
			cout << "Forbidding subtractions\n";
		else
			cerr << "Allowing subtractions\n";

		double qsum = 0.0;
		for (int i = 1; i <= K; i++) {
			Q0[i] = pu[i-1];
			qsum += Q0[i];
		}
		for (int i = 1; i <= K; i++) {
			if (qsum != 1.0)
				Q0[i] /= qsum;
			lgQ0[i] = log(Q0[i]);
		}

		// calculate recommended permutation for convex
		for (int i = 1; i <= K; i++)
			perm[i] = i;
		std::sort( perm + 1, perm + 1 + K, perm_cmp( Q0 ) );
		//qsort(perm + 1, K, sizeof(int), compare_forperm);

		for (int i = 1; i <= K; i++) {
			Q[i] = Q0[perm[i]];
			lgQ[i] = lgQ0[perm[i]]; // = log(Q[i])
		}
		double qcomp = 1.0;
		lgQcomp[0] = 0.0; // boundary condition for reverse Q sum
		for (int i = 1; i <= K - 1; i++) {
			Qcomp[i] = (qcomp -= Q[i]);
			Qcompmin[i] = d_min(Q, i + 1, K); // all Q is now assigned
			lgQcomp[i] = log(Qcomp[i]);
			lgQcompmin[i] = log(Qcompmin[i]); // all Q is now assigned
		}

		lg = cumlg = 0;
		biggest_N = 0;
	}

	~pval_pvalue_calculator() {
		delete[] Q0;
		delete[] lgQ0;
		delete[] C0;
		delete[] n;
		delete[] perm;
		delete[] Q;
		delete[] lgQ;
		delete[] lgQcomp;
		delete[] lgQcompmin;
		delete[] Qcomp;
		delete[] Qcompmin;

	#ifdef LOOKUP
		delete[] logaddf;
		delete[] logsubf;
	#endif

		delete[] lg;
		delete[] cumlg;
	}

	double
	operator()( const p_value_test_case & args ) {
		// update N and resize our log[i] arrays if needed
		set_N( args.N );

		// calculate the p-value
		double log_pop = 1.;
		BOOST_FOREACH( double IC, args.ICs ) {
			log_pop += calc_col_pvalue( IC );
		}
		return use_qfast ? log_qfast( args.W(), log_pop ) : log_pop;
	}

	void
	set_N( size_t new_N ) {
		N = new_N;
		if( N > biggest_N ) {
			delete [] lg;
			delete [] cumlg;
			biggest_N = N;

			lg = new double[N + 1]; // lg[i]    = log(i)
			cumlg = new double[N + 1]; // cumlg[i] = log(i!)
			lg[0] = cumlg[0] = 0.0; // cumlg[0] inits the sum
			// both used as boundary conditions for recursions so that 0*log(0) = 0
			for (int i = 1; i <= N; i++)
				cumlg[i] = (cumlg[i - 1] + (lg[i] = log(i)));
		}
	}

	double
	calc_col_pvalue( double IC ) {
		eps = IC * N;
		KL_thresh = eps + N * lg[N];

		reccalls = logops = discards = 0; // updated only in cases 1,2,3
		underflow = 0; // updated only in 3
		clock_t start, finish;
		time_t start_sec, finish_sec; // in case compute time too long for clock()
		int df;

		switch (alg) {
		case 0:
			df = K - 1;
			start_sec = time(NULL);
			start = clock();
			pval = chi2pval(df, 2 * eps);
			lg_pval = log(pval);
			finish = clock();
			finish_sec = time(NULL);
			break;

		case 1:
		case 2:
		case 3:
			lg_pval = -HUGE_VAL; // ==log(0)
			start_sec = time(NULL);
			start = clock();
			if (alg == 1)
				exact(1, 0, 0.0, cumlg[N]);
			else if (alg == 2)
				expand(1, 0, 0.0, cumlg[N]);
			else
				convex(1, 0, 0.0, cumlg[N]);
			//pval = exp(lg_pval);
			finish = clock();
			finish_sec = time(NULL);
			break;

		case 4:
			throw std::logic_error("Simulation not implemented.");
			start_sec = time(NULL);
			start = clock();
			pval = g2pvalsim(simsec);
			lg_pval = log(pval);
			finish = clock();
			finish_sec = time(NULL);
			break;

		default:
			cout << "Got unknown algorithm argument.\n";
			exit(1);
		}

		runtime = double(finish - start) / CLOCKS_PER_SEC;

		return lg_pval;

	}

	struct perm_cmp {
		double * Q0;
		perm_cmp( double * Q0 ) : Q0( Q0 ) { }
		bool operator()( int i, int j ) {
			return Q0[ i ] < Q0[ j ];
		}
	};

	void exact(int i, int n_sum, double KL_sum, double Qtau_sum) {// we are guaranteed that lg[0]=cumlg[0]=0.0

		reccalls++;
		for (n[i] = N - n_sum; n[i] >= 0; n[i]--) // never descends to level K
		{
			double KL_term = n[i] * (lg[n[i]] - lgQ[i]);
			double Qtau_add = n[i] * lgQ[i] - cumlg[n[i]];

			if (i < K - 1)
				exact(i + 1, n_sum + n[i], KL_sum + KL_term, Qtau_sum + Qtau_add);
			else { // its faster not to recurse to level K
				n[K] = N - n_sum - n[i];
				KL_term += n[K] * (lg[n[K]] - lgQ[K]);
				Qtau_add += n[K] * lgQ[K] - cumlg[n[K]];
				if (KL_sum + KL_term >= KL_thresh) {
					logops++;
					add_lgp(Qtau_sum + Qtau_add);
				} else
					discards++;
			}

		}
	}

	void expand(int i, int n_sum, double KL_sum, double Qtau_sum) { // we are guaranteed that lg[0]=cumlg[0]=0.0

		reccalls++;
		for (n[i] = N - n_sum; n[i] >= 0; n[i]--) // never descends to level K
		{
			double KL_term = n[i] * (lg[n[i]] - lgQ[i]);
			int n_comp = N - n_sum - n[i];
			double Qtau_add = n[i] * lgQ[i] - cumlg[n[i]];
			double KL_minlf = n_comp * (lg[n_comp] - lgQcomp[i]);

			if (KL_sum + KL_term + KL_minlf >= KL_thresh) {
				logops++;
				add_lgp(Qtau_sum + Qtau_add + n_comp * lgQcomp[i] - cumlg[n_comp]);
			} else if (i < K - 1) {
				double KL_maxlf = n_comp * (lg[n_comp] - lgQcompmin[i]);
				if (KL_sum + KL_term + KL_maxlf >= KL_thresh)
					expand(i + 1, n_sum + n[i], KL_sum + KL_term, Qtau_sum
							+ Qtau_add);
				else
					// all subtrees below thresh
					discards++;
			} else
				// lgQcomp[K-1]==lgQcompmin[K-1]=lgQ[K] thus minlf==maxlf
				discards++;

		}
	}

	void convex(int i, int n_sum, double KL_sum, double Qtau_sum) {// we are guaranteed that lg[0]=cumlg[0]=0.0

		reccalls++;
		int n_comp = N - n_sum;
		int lmin, rmin, lmax, rmax, a, b, c, d, dummy, mmin, mmax;

		mmin = int(Q[i] * (n_comp) / (Q[i] + Qcomp[i]));
		a = 0;
		dummy = mmin;
		lmin = binsearch_descend(a, dummy, KL_thresh - KL_sum, n_comp, i, lgQcomp);
		dummy = mmin + 1;
		d = n_comp;
		rmin = binsearch_ascend(dummy, d, KL_thresh - KL_sum, n_comp, i, lgQcomp);

		mmax = int(Q[i] * (n_comp) / (Q[i] + Qcompmin[i]));
		dummy = 0;
		b = mmax;
		lmax = binsearch_descend(dummy, b, KL_thresh - KL_sum, n_comp, i,
				lgQcompmin);
		c = mmax + 1;
		dummy = n_comp;
		rmax
				= binsearch_ascend(c, dummy, KL_thresh - KL_sum, n_comp, i,
						lgQcompmin);

		if (lmin == -1 || (lmin == 1 && rmin == 1))
			a = -1;
		else if (lmin == 1)
			a = mmin;

		if (rmin == -1)
			d = n_comp + 1;
		else if (rmin == 1)
			d = (lmin == 1 ? 0 : mmin + 1);

		if (lmax == -1)
			b = 0;
		else if (lmax == 1)
			b = (rmax == 1 ? n_comp + 1 : mmax + 1);

		if (rmax == -1)
			c = n_comp;
		else if (rmax == 1)
			c = (lmax == 1 ? -1 : mmax);

		if (FORBID_SUBTRACTIONS || n_comp + 3 <= 2 * (d - a))
		// if (logops: |add| = a+2+n_comp-d <= |sub| = d-a-1)
		{
			logops += a + 2 + n_comp - d;
			for (n[i] = n_comp; n[i] >= 0; n[i]--) // never descends to level K
			{
				if (n[i] <= a || n[i] >= d) // (add)
					add_lgp(Qtau_sum + n[i] * lgQ[i] - cumlg[n[i]]
							+ (n_comp - n[i]) * lgQcomp[i] - cumlg[n_comp - n[i]]);
				else if (i < K - 1 && (n[i] < b || n[i] > c)) // (descend)
					// lgQcomp[K-1]==lgQcompmin[K-1]=lgQ[K] thus not add => discard
					convex(i + 1, n_sum + n[i],
							KL_sum + n[i] * (lg[n[i]] - lgQ[i]), Qtau_sum + n[i]
									* lgQ[i] - cumlg[n[i]]);
				else
					// (discard)
					discards++;
			}
		} else // faster (but less accurate) to subtract
		{ // as written - do NOT mix with approx pval. overshoots, only LATER mends!
			logops += d - a - 1;
			add_lgp(Qtau_sum + n_comp * lgQcomp[i - 1] - cumlg[n_comp]);
			for (n[i] = d - 1; n[i] >= a + 1; n[i]--) // never descends to level K
			{ // (loop over all not add)
				sub_lgp(Qtau_sum + n[i] * lgQ[i] - cumlg[n[i]] + (n_comp - n[i])
						* lgQcomp[i] - cumlg[n_comp - n[i]]);
				if (i < K - 1 && (n[i] < b || n[i] > c)) // (descend)
					// lgQcomp[K-1]==lgQcompmin[K-1]=lgQ[K] thus not add => discard
					convex(i + 1, n_sum + n[i],
							KL_sum + n[i] * (lg[n[i]] - lgQ[i]), Qtau_sum + n[i]
									* lgQ[i] - cumlg[n[i]]);
				else
					// (discard)
					discards++;
			}
		}
	}

	void convex_check(int i, int n_sum, double KL_sum, double Qtau_sum) {// we are guaranteed that lg[0]=cumlg[0]=0.0
		// we run exact and check every node for the convex conditions
		// debugging routine.

		int n_comp = N - n_sum;

		if (i < K) {
			int lmin, rmin, lmax, rmax, a, b, c, d, dummy, mmin, mmax;
			mmin = int(Q[i] * (n_comp) / (Q[i] + Qcomp[i]));
			a = 0;
			dummy = mmin;
			lmin = binsearch_descend(a, dummy, KL_thresh - KL_sum, n_comp, i,
					lgQcomp);
			dummy = mmin + 1;
			d = n_comp;
			rmin = binsearch_ascend(dummy, d, KL_thresh - KL_sum, n_comp, i,
					lgQcomp);

			mmax = int(Q[i] * (n_comp) / (Q[i] + Qcompmin[i]));
			dummy = 0;
			b = mmax;
			lmax = binsearch_descend(dummy, b, KL_thresh - KL_sum, n_comp, i,
					lgQcompmin);
			c = mmax + 1;
			dummy = n_comp;
			rmax = binsearch_ascend(c, dummy, KL_thresh - KL_sum, n_comp, i,
					lgQcompmin);

			if (lmin == -1 || (lmin == 1 && rmin == 1))
				a = -1;
			else if (lmin == 1)
				a = mmin;

			if (rmin == -1)
				d = n_comp + 1;
			else if (rmin == 1)
				d = (lmin == 1 ? 0 : mmin + 1);

			if (lmax == -1)
				b = 0;
			else if (lmax == 1)
				b = (rmax == 1 ? n_comp + 1 : mmax + 1);

			if (rmax == -1)
				c = n_comp;
			else if (rmax == 1)
				c = (lmax == 1 ? -1 : mmax);

			if (a >= b || c >= d)
				cerr << "(" << lmin << "," << rmin << "|" << lmax << "," << rmax
						<< ")";
		}

		for (n[i] = n_comp; n[i] >= 0; n[i]--) // evals all at level K
		{
			reccalls++;
			double KL_term = n[i] * (lg[n[i]] - lgQ[i]);
			double Qtau_add = n[i] * lgQ[i] - cumlg[n[i]];

			if (i < K)
				convex(i + 1, n_sum + n[i], KL_sum + KL_term, Qtau_sum + Qtau_add);
			else if (KL_sum + KL_term >= KL_thresh) {
				logops++;
				add_lgp(Qtau_sum + Qtau_add);
			} else
				discards++;

			if (i == K) // single allowed level K assignment
				break;
		}
	}

	int binsearch_descend(int &l, int &r, double c, int n_comp, int i,
			double *lgq_star) { // input : r>=l and the implicit F is monotonically decreasing
		// output: -1  F(all)<c
		//          0  F(l)>=c>f(r)
		//         +1  F(all)>=c

		if (F(l) < c)
			return (-1);
		if (F(r) >= c)
			return (+1);
		while (r - l > 1) {
			int m = (l + r) / 2;
			if (F(m) < c)
				r = m;
			else
				l = m;
		}
		return (0);
	}

	int binsearch_ascend(int &l, int &r, double c, int n_comp, int i,
			double *lgq_star) { // input : r>=l and the implicit F is monotonically increasing
		// output: -1  F(all)<c
		//          0  F(l)<c<=f(r)
		//         +1  F(all)>=c
		// (couldn't use descending with -F because the equality changes place)

		if (F(r) < c)
			return (-1);
		if (F(l) >= c)
			return (+1);
		while (r - l > 1) {
			int m = (l + r) / 2;
			if (F(m) < c)
				l = m;
			else
				r = m;
		}
		return (0);
	}

	#undef F

	// inline double addlgs (double lgx, double lgy)
	// {
	//   // input:  lg(x), lg(y)
	//   // output: lg(x+y)
	//   // note: log(0) == -HUGE_VAL so use -HUGE_VAL to init sums!!!

	//   if (lgx == -HUGE_VAL)
	//     return(lgy);
	//   if (lgy == -HUGE_VAL)
	//     return(lgx);

	//   if (lgx < lgy)
	//     return(lgy+log(1.0+exp(lgx-lgy)));

	//   return(lgx+log(1.0+exp(lgy-lgx)));
	// }


	inline void add_lgp(double lgy) {
		// global lg_pval=log(pval) + input lgy=log(y)  => lg_pval=log(pval+y)
		// note: log(0) == -HUGE_VAL so use -HUGE_VAL to init sums!!!

	#ifdef LOOKUP
		if (lg_pval == -HUGE_VAL)
		lg_pval = lgy;
		else if (lgy != -HUGE_VAL)
		{
			if (lg_pval < lgy)
			{
				double dummy = lg_pval;
				lg_pval = lgy;
				lgy = dummy;
			}
			int i = int((lg_pval-lgy)*G_LOGADD);
			if (i<d_logadd) // linear interpolation
			lg_pval += ((i+1-(lg_pval-lgy)*G_LOGADD)*logaddf[i] +
					((lg_pval-lgy)*G_LOGADD-i)*logaddf[i+1]);
		}
		return;
	#else
		if (lg_pval == -HUGE_VAL)
			lg_pval = lgy;
		else if (lgy != -HUGE_VAL) {
			if (lg_pval < lgy)
				lg_pval = lgy + log(1.0 + exp(lg_pval - lgy));
			else
				lg_pval = lg_pval + log(1.0 + exp(lgy - lg_pval));
		}
		return;
	#endif
	}

	// inline double sublgs (double lgx, double lgy)
	// {
	//   // input:  lg(x) > lg(y) >= log(0) == -HUGE_VAL
	//   // output: lg(x-y)

	//   if (lgx <= lgy)
	//     {
	//       underflow++;
	//       return(-HUGE_VAL);
	//     }

	//   if (lgy == -HUGE_VAL)
	//     return(lgx);

	//   return(lgx+log(1.0-exp(lgy-lgx)));
	// }

	inline void sub_lgp(double lgy) {
		// global lg_pval=log(pval) + input lgy=log(y)  => lg_pval=log(pval-y)
		//   or log(0)=-HUGE_VAL if pval<=y

	#ifdef LOOKUP
		if (lg_pval <= lgy)
		{
			underflow++;
			lg_pval = -HUGE_VAL;
		}
		else if (lgy != -HUGE_VAL)
		if (lg_pval-lgy < .02) // too steep to interpolate on
		lg_pval = lg_pval+log(1.0-exp(lgy-lg_pval));
		else
		{
			/*{
			 int i = int(-lgy);
			 if (i<=E_BINS)
			 e_bin[i]++;
			 }*/
			int i = int((lg_pval-lgy)*G_LOGSUB);
			if (i<d_logsub) // linear interpolation
			lg_pval += ((i+1-(lg_pval-lgy)*G_LOGSUB)*logsubf[i] +
					((lg_pval-lgy)*G_LOGSUB-i)*logsubf[i+1]);
		}

		return;
	#else
		if (lg_pval <= lgy) {
			underflow++;
			lg_pval = -HUGE_VAL;
		} else if (lgy != -HUGE_VAL)
			lg_pval = lg_pval + log(1.0 - exp(lgy - lg_pval));

		return;
	#endif
	}

	double g2pvalsim(double secs) {
		int rounds;
		logops = 0;

		double *cumQ = new double[K + 1];
		cumQ[0] = 0.0;
		for (int i = 1; i <= K; i++)
			cumQ[i] = cumQ[i - 1] + Q0[i] * RAND_MAX;

		time_t start = clock();
		for (rounds = 1; double(clock() - start) / CLOCKS_PER_SEC < secs; rounds++) {
			for (int i = 1; i <= K; i++)
				n[i] = 0;

			for (int j = 1; j <= N; j++) //generate a sample of size N
			{
				int num = rand();
				for (int i = 1; i <= K; i++)
					if (num <= cumQ[i]) {
						n[i]++;
						break;
					}
			}

			double halfg2 = -N * lg[N];
			for (int i = 1; i <= K; i++) {
				if (n[i] > 0)
					halfg2 += n[i] * (lg[n[i]] - lgQ0[i]);
			}
			if (halfg2 >= eps)
				logops++;
		}

		delete[] cumQ;
		discards = rounds - logops;
		return (double(logops) / rounds);
	}

};


pvalue_calculator::ptr
create_PVAL_pvalue_calculator(
	double * pu,
	int alg,
	bool use_qfast
) {
	return pvalue_calculator::ptr( new pval_pvalue_calculator( alg, pu, use_qfast ) );
}





//-----------begin NR code ---------------------

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}

double gammln(double xx) {
	double x, y, tmp, ser;
	static double cof[6] = { 76.18009172947146, -86.50532032941677,
			24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
			-0.5395239384953e-5 };
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j <= 5; j++)
		ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005 * ser / x);
}

#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln) {
	double gammln(double xx);
	void nrerror(char error_text[]);
	int n;
	double sum, del, ap;

	*gln = gammln(a);
	if (x <= 0.0) {
		if (x < 0.0)
			nrerror( const_cast< char * >( "x less than 0 in routine gser" ));
		*gamser = 0.0;
		return;
	} else {
		ap = a;
		del = sum = 1.0 / a;
		for (n = 1; n <= ITMAX; n++) {
			++ap;
			del *= x / ap;
			sum += del;
			if (fabs(del) < fabs(sum) * EPS) {
				*gamser = sum * exp(-x + a * log(x) - (*gln));
				return;
			}
		}
		nrerror( const_cast< char * >( "a too large, ITMAX too small in routine gser" ) );

		return;
	}
}

#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln) {
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an, b, c, d, del, h;

	*gln = gammln(a);
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;
	for (i = 1; i <= ITMAX; i++) {
		an = -i * (i - a);
		b += 2.0;
		d = an * d + b;
		if (fabs(d) < FPMIN)
			d = FPMIN;
		c = b + an / c;
		if (fabs(c) < FPMIN)
			c = FPMIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.0) < EPS)
			break;
	}
	if (i > ITMAX)
		nrerror( const_cast< char * >( "a too large, ITMAX too small in gcf" ) );
	*gammcf = exp(-x + a * log(x) - (*gln)) * h;
}
#undef ITMAX
#undef EPS
#undef FPMIN

double gammq(double a, double x) {
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser, gammcf, gln;

	if (x < 0.0 || a <= 0.0)
		nrerror( const_cast< char * >( "Invalid arguments in routine gammq" ) );
	if (x < (a + 1.0)) {
		gser(&gamser, a, x, &gln);
		return 1.0 - gamser;
	} else {
		gcf(&gammcf, a, x, &gln);
		return gammcf;
	}
}

//-----------end NR code ---------------------



