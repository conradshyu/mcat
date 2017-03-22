/*
 * strain.h
 *
 * Written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * Center for the Study of Biological Complexity (CSBC)
 * Department of Microbiology and Immunology
 * Medical College of Virginia
 * Virginia Commonwealth University
 * Richmond, VA 23298
 *
 * revised on April 15, 2013
*/

#ifndef _STRAIN_H
#define _STRAIN_H

#include <table.h>
#include <pivot.h>

#include <map>
#include <vector>
#include <string>

class Strain
{
public:
    Strain( const std::map<unsigned int, stTABLE>& );
    ~Strain();

    bool Run( const std::string& );
    const std::map<unsigned int, double>& GetIndex() const;

private:
    std::map<unsigned int, double> mIndex;
    std::map<unsigned int, stPIVOT> mAssign;
    std::map<unsigned int, std::string> mTaxon;
    std::map<unsigned int, unsigned int> mBlock;

    bool Assign( const std::string& );
    bool Output( const std::string& );

    double Weight( const std::vector<unsigned int>&, const double ) const;
    double Shannon( const std::vector<unsigned int>&, const double = 1.0 ) const;
};  // end of class definition

#endif  // _ASSIGN_H
