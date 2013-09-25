#ifndef _LIST_H_
#define _LIST_H_

//
// Class for storing the values of the pmf
// in a list (for Hirji's algorithm). 
// The values can also be stored in a array (for e.g. array.h,
// not included in this package) or any other data structure 
// which has the same public interface as this class.
//
// elem_t = datatype for the array elements
//
template< class elem_t >
    class list_t {

    private:

        //
        // The "list" is actually maintained as an array of indices for
        // non-zero values of the list (index) and an array of indexed
        // values of the list (data).  The improvement this data structure
        // provides over arrays is that it allows for fast iteration over
        // the non-zero values of the array.
        //
        // data[0..size-1] contains the data.
        // cur_i is an index for iterating over the array.
        // index[0..max_index-1] contain the indices for the non-zero values.
        //
        elem_t* data;
        int size;
        int cur_i;
        int* index;
        int max_index;

    public:

        list_t() {
            size = 0;
        }

        void
        init(
            int given_size ) {

            size = given_size;
            data = new elem_t[ size ];
            std::fill(
                data, data + size, 0. );
            max_index = 0;
            index = new int[ size + 1 ];
        }

        list_t(
            int given_size ) {

            init(
                given_size );
        }

        inline void
        set(
            int i, elem_t value ) {

            if( data[ i ] == 0 )
                index[ max_index++ ] = i;
            data[ i ] = value;
        }

        inline elem_t&
        operator[](
            int i ) {
            return data[ i ];
        }

        //
        // Sets the value at the current index to 0
        //
        inline void
        remove_cur() {

            data[ index[ cur_i ] ] = 0;
            index[ cur_i ] = index[ --max_index ];
            cur_i--;
        }

        //
        // Functions for iterating over the non-zero values of the array
        // using the variable cur_i.
        //
        // start: Initializes cur_i and returns the first index.
        // not_last: Returns 1 if cur_i is not at the last index.
        // next: Returns the next index corresponding to a non-zero value.
        //

        //
        // minI is a dummy parameter so the signature matches that
        // of array.h
        //
        inline int
        start(
            int minI = 0 ) {
            cur_i = 0;
            return index[ cur_i ];
        }

        //
        // maxI is a dummy parameter so the signature matches that
        // of array.h
        //
        inline int
        not_last(
            int maxI = 0 ) {
            return cur_i < max_index;
        }

        //
        // maxI is a dummy parameter so the signature matches that
        // of array.h
        //
        inline int
        next(
            int maxI = 0 ) {
            return index[ ++cur_i ];
        }

        inline void
        destructive_copy(
            list_t< elem_t >& a ) {
            delete[] ( data );
            delete[] ( index );
            *this = a;
            a.data = 0;
            a.index = 0;
        }

        ~list_t() {
            //delete[](data); delete[](index); data = 0; index = 0;
        }

    };

#endif
