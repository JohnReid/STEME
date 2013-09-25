/**
 * Copyright John Reid 2011
 *
 * \file
 * \brief Test memory usage of std::vector.
 */

#include <boost/test/utils/wrap_stringstream.hpp>

#include <vector>
#include <iostream>
#include <stdexcept>

#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>

#define MAKE_STRING( x ) ( boost::wrap_stringstream().ref() << x ).str()

using namespace std;

template< typename T, class Alloc >
void
shrink_to_empty( vector< T, Alloc > & v )
{
    // swap with empty vector
    vector< T, Alloc >( v.begin(), v.end(), v.get_allocator() ).swap( v );
}

void
run_command( const std::string & command ) {
    if( system( command.c_str() ) ) {
        throw std::logic_error( MAKE_STRING( "Could not run command: " << command ) );
    }
}

int
main( int argc, char * argv[] ) {
    std::vector< int > v;
    char buf[2048];
    const size_t size = 1e8;
    const std::string show_process_info_cmd = MAKE_STRING( "ps -Fp " << getpid() );
    run_command( show_process_info_cmd );

    v.reserve( size );
    cout << "Reserved " << size << " elements.\n";
    run_command( show_process_info_cmd );

    v.resize( size );
    cout << "Resized " << size << " elements.\n";
    run_command( show_process_info_cmd );

    v.clear();
    cout << "Cleared " << size << " elements.\n";
    run_command( show_process_info_cmd );

    v.reserve( 0 );
    cout << "Reserved 0 elements.\n";
    run_command( show_process_info_cmd );

    shrink_to_empty( v );
    cout << "Shrunk to empty.\n";
    run_command( show_process_info_cmd );

    return 0;
}
