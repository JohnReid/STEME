#include "utils.h"
#include "llr_score.h"

//
// Initializes all the data members
// Q = Lattice size
//
llr_score_t::llr_score_t(
    int given_N, int given_K, double* given_pu, int given_Q ) {

    N = given_N;
    K = given_K;
    pu = given_pu;
    Q = given_Q;

    double* log_pu = new double[ K ];
    for( int i = 0; i < K; i++ )
        log_pu[ i ] = log(
            pu[ i ] );

    double* logn = new double[ N + 1 ];
    logn[ 0 ] = LOGZERO;
    for( int i = 1; i < N + 1; i++ )
        logn[ i ] = log(
            ( double )i );

    int imin;
    min_array(
        pu, K, imin );
    max_ent = -N * log_pu[ imin ];
    step = max_ent / ( Q - 1 );
    min_t_ent = -N * ( 1 - pu[ imin ] ) / exp(
        1.0 );

    //
    // Ip[k][n] = round(n*log(n/(N*pu[k]))/step)
    //
    Ip = new int*[ K ];
    I = new double*[ K ];
    for( int i = 0; i < K; i++ ){

        Ip[ i ] = new int[ N + 1 ];
        Ip[ i ][ 0 ] = 0;
        I[ i ] = new double[ N + 1 ];
        I[ i ][ 0 ] = 0.0;
        for( int j = 1; j < N + 1; j++ ){
            I[ i ][ j ] = j * ( logn[ j ] - logn[ N ] - log_pu[ i ] );
            Ip[ i ][ j ] = NINT(I[i][j]/step);
        }
    }

    double min_ent_N[ N + 1 ];
    for( int n = 0; n < N + 1; n++ )
        min_ent_N[ n ] = I[ 0 ][ n ];

    for( int k = 1; k < K; k++ )
        for( int n = N; n >= 0; n-- )
            for( int j = n; j >= 0; j-- )
                min_ent_N[ n ] = MIN(min_ent_N[n], min_ent_N[j]+Ip[k][n-j]);

    min_ent = min_ent_N[ N ];

    delete[] ( log_pu );
    delete[] ( logn );
}

llr_score_t::~llr_score_t() {

    for( int i = 0; i < K; i++ )
        delete[] ( Ip[ i ] );
    delete[] ( Ip );
}
