#include "basic_pruning.h"
#include "utils.h"
#include <iostream>
using namespace std;

//
// s = Score for which pruning is carried out
//
template< class score_t, class prob_t, class multi_pmf_t >
    basic_pruning_t< score_t, prob_t, multi_pmf_t >::basic_pruning_t(
        score_t* given_score, prob_t* given_prob, multi_pmf_t* given_multi_pmf,
        int offset, double s ) {

        score = given_score;
        multi_pmf = given_multi_pmf;
        int N = score->get_N();
        int K = score->get_K();
        kn2 = 1 / ( ( double )K * N * ( N - 1 ) / 2.0 );
        Qavg = 0;

        //
        // Initialization (See the .h file for the definition of min_firstK and max_firstK)
        //
        min_firstK_I = new int*[ K ];
        max_firstK_I = new int*[ K ];
        min_firstK_I[ 0 ] = new int[ N + 1 ];
        max_firstK_I[ 0 ] = new int[ N + 1 ];
        for( int n = 0; n < N + 1; n++ )
            min_firstK_I[ 0 ][ n ] = max_firstK_I[ 0 ][ n ] = score->get_Ip(
                0, n ) - offset;

        //
        // Recursively computing the bounds
        //

        //
        // Keeps track of the minimum index seen so far.
        //
        int min_I = 0;

        for( int k = 1; k < K; k++ ){

            min_firstK_I[ k ] = new int[ N + 1 ];
            max_firstK_I[ k ] = new int[ N + 1 ];
            for( int n = N; n >= 0; n-- ){

                min_firstK_I[ k ][ n ] = min_firstK_I[ k - 1 ][ n ]
                                + score->get_Ip(
                                    k, 0 );
                max_firstK_I[ k ][ n ] = max_firstK_I[ k - 1 ][ n ]
                                + score->get_Ip(
                                    k, 0 );

                //
                // i chooses the number of objects in the kth bin.
                //
                for( int i = 1; i <= n; i++ ){
                    min_firstK_I[ k ][ n ] = MIN(min_firstK_I[k][n],
                                    min_firstK_I[k-1][n-i]+score->get_Ip(k, i));
                    max_firstK_I[ k ][ n ] = MAX(max_firstK_I[k][n],
                                    max_firstK_I[k-1][n-i]+score->get_Ip(k, i));
                }

                if( min_firstK_I[ k ][ n ] < min_I )
                    min_I = min_firstK_I[ k ][ n ];
            }
        }

#if DEBUG
        if(min_I < 0)
        cout << "Negative indices: " << min_I << endl;
#endif

    }
