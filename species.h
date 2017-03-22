/*
 * species.h
 *
 * Written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * Center for the Study of Biological Complexity (CSBC)
 * Department of Microbiology and Immunology
 * Medical College of Virginia
 * Virginia Commonwealth University
 * Richmond, VA 23298
 *
 * revised on April 16, 2013
*/

#ifndef _SPECIES_H
#define _SPECIES_H

#include <table.h>
#include <pivot.h>

#include <map>
#include <string>

class Species
{
public:
    Species(
        const std::map<unsigned int, stTABLE>&,     // translate table
        const std::map<unsigned int, double>& );    // weighted shannon index
    ~Species();

    bool Run( const std::string& );

private:
    std::map<std::string, double> mIndex;
    std::map<std::string, stPIVOT> mAssign;
    std::map<unsigned int, std::string> mTaxon;

    bool Output( const std::string& );
    bool Profile( const std::string& );
    bool Assign( const std::string& );
    bool Assign( const std::string&, const stPIVOT& );
};  // end of class definition

#endif  // _SPECIES_H
