#ifndef _MULTI_PMF_H_
#define _MULTI_PMF_H_

//
// Class for an array of pmfs of type pmf_t.
// This data structure is maintained as
// part of Hirji's DP algorithm. Strictly speaking,
// the nth "pmf" typically needs to be multiplied by a function
// of n in order to be a pmf; for e.g. for the multinomial
// distribution this function is f(n)=n!
//
template< class pmf_t >
    class multi_pmf_t {

    private:

        //
        // N = Number of pmfs
        // max_range = Size of pmfs
        //
        int N;
        int max_range;

    public:

        //
        // pmf[n] = The nth pmf
        //
        pmf_t* pmf;

        multi_pmf_t(
            int given_N, int given_max_range ) {
            N = given_N;
            max_range = given_max_range;
            pmf = new pmf_t[ N ];
            for( int i = 0; i < N; i++ )
                pmf[ i ].init(
                    max_range );
        }

        inline int
        get_max_range() {
            return max_range;
        }

        //
        // Copies the pmf from c while setting its pointer to 0.
        //
        inline void
        destructive_copy(
            multi_pmf_t< pmf_t >& c ) {
            delete[] ( pmf );
            *this = c;
            c.pmf = 0;
        }

        ~multi_pmf_t() {
            delete[] ( pmf );
            pmf = 0;
        }

    };

#endif
