/*
 * assign.cpp
 *
 * Written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * Center for the Study of Biological Complexity (CSBC)
 * Department of Microbiology and immunology
 * Medical College of Virginia
 * Virginia Commonwealth University
 * Richmond, VA 23298
 *
 * revised on April 16, 2013
*/

#include <table.h>
#include <strain.h>
#include <species.h>

#include <fstream>
#include <iostream>

/*
 * main driver procedure
*/
int main( int argc, char* argv[] )
{
    const unsigned int nMaxBUFFER = 2048;

    if ( argc < 3 )
    {
        return( 1 );
    }   // check the number of parameters

    stTABLE a;
    char buffer[ nMaxBUFFER ];
    std::map<unsigned int, stTABLE> table;

    std::ifstream ifs( argv[ 1 ], std::ios::in );

    if ( ifs.fail() )
    {
        return( 0 );
    }   // check the state of stream

    ifs.getline( buffer, nMaxBUFFER );      // skip the header
    std::cout << "loading translation table ..." << std::flush;

    while ( ifs.getline( buffer, nMaxBUFFER ) )
    {
        a = buffer; table[ a.gid ] = a;
    }   // load the table first

    ifs.close(); std::cout << " completed" << std::endl;

    std::cout << "processing file: " << argv[ 2 ] << std::endl;
    std::cout << "strain level assignment ..." << std::flush;
    Strain p( table );                  // strain level assignment
    p.Run( argv[ 2 ] );
    std::cout << " completed" << std::endl;

    std::cout << "species level assignment ..." << std::flush;
    Species q( table, p.GetIndex() );   // species level assignment
    q.Run( argv[ 2 ] );
    std::cout << " completed" << std::endl;

    return( 0 );
}   // end of main()
