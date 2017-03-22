/*
 * table.h
 *
 * Written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * Center for the Study of Biological Complexity
 * Department of Microbiology and Immunology
 * Medical College of Virgina
 * Virginia Commonwealth University
 * Richmond, VA 23298
 *
 * load the translation table for further processing
 * smart container implementation
 *
 * revised on January 8, 2013
 * revised on January 9, 2013
 * revised on February 5, 2013
 * revised on February 26, 2013
 * revised on April 10, 2013
*/

#ifndef _TABLE_H   // only load it once
#define _TABLE_H

// c++ specific headers
#include <vector>
#include <string>
#include <cstring>
#include <boost/algorithm/string.hpp>

struct stTABLE
{
    stTABLE( const stTABLE& _t )
    {
        *this = _t;
    }   // copy constructor

    stTABLE( const std::string& _s )
    {
        *this = _s;
    }   // assignment constructor

    stTABLE()
    {
        strain.clear(); species.clear();
    }   // default constructor

    ~stTABLE()
    {
        strain.clear(); species.clear();
    }  // default destructor; environmentally conscientious

    const stTABLE& operator=( const stTABLE& _t )
    {
        gid = _t.gid, tid = _t.tid;     // ncbi gid and tid
        start = _t.start, end = _t.end; // start and end of histogram bin
        size = _t.size;                 // size of genome
        strain = _t.strain;             // strain name
        species = _t.species;           // species name

        return( *this );
    }   // end of operator overloading

    const stTABLE& operator=( const std::string& _s )
    {
        const char* szTOKEN = ",\n";
        std::vector<std::string> field;

        boost::algorithm::split(                // splite the entire string
            field, _s, boost::algorithm::is_any_of( szTOKEN ) );

        gid = static_cast<unsigned int>( ::atoi( field[ 0 ].c_str() ) );    // ncbi gid
        tid = static_cast<unsigned int>( ::atoi( field[ 1 ].c_str() ) );    // ncbi tid
        size = static_cast<unsigned int>( ::atoi( field[ 2 ].c_str() ) );   // size of genome
        start = static_cast<unsigned int>( ::atoi( field[ 3 ].c_str() ) );  // start of bin
        end = static_cast<unsigned int>( ::atoi( field[ 4 ].c_str() ) );    // end of bin
        strain = field[ 5 ]; species = field[ 6 ];
/*
        // keep the quotation marks around strings?
        boost::algorithm::trim_if( strain, boost::algorithm::is_any_of( "\"" ) );
        boost::algorithm::trim_if( species, boost::algorithm::is_any_of( "\"" ) );
*/
        return( *this );
    }   // end of operator overloading

    unsigned int gid;       // ncbi gid
    unsigned int tid;       // ncbi tid
    unsigned int size;      // size of genome
    unsigned int start;     // start of histogram bin
    unsigned int end;       // end of histogram bin
    std::string strain;     // strain name
    std::string species;    // species name
};  // smart container

#endif  // _TABLE_H
