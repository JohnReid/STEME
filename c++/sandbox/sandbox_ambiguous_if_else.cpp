/**
 * Copyright John Reid 2011
 */

#include <iostream>

using namespace std;

void
test_ambiguous( bool a, bool b, bool c ) {
	// to test how C interprets ambiguous else
	cout << "AMBIGUOUS:     " << a << " " << b << " " << c << ": ";
	if( a )
	{
		cout << "1";
	}
	else if( b )
	if( c )
	cout << "2";
	else
	{
		cout << "3";
	}
	cout << "\n";
}


void
test_disambiguated( bool a, bool b, bool c ) {
	// to test how C interprets ambiguous else
	cout << "DISAMBIGUATED: " << a << " " << b << " " << c << ": ";
	if( a )
	{
		cout << "1";
	} else if( b ) {
		if( c ) {
			cout << "2";
		} else {
			cout << "3";
		}
	}
	cout << "\n";
}

void
test( bool a, bool b, bool c ) {
	test_ambiguous( a, b, c );
	test_disambiguated( a, b, c );
}

int
main( int argc, char * argv[] ) {

	test( true, true, true );
	test( true, true, false );
	test( true, false, true );
	test( true, false, false );
	test( false, true, true );
	test( false, true, false );
	test( false, false, true );
	test( false, false, false );

	return 0;
}
