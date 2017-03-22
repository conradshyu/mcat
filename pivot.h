/*
 * pivot.h
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

#ifndef _PIVOT_H
#define _PIVOT_H

#include <vector>

struct stPIVOT
{
    stPIVOT()
    {
        site.clear();
    }   // default constructor

    stPIVOT( const stPIVOT& _t )
    {
        *this = _t;
    }   // copy constructor

    ~stPIVOT()
    {
        site.clear();
    }   // default destructor; environmentally conscientious

    const stPIVOT& operator=( const stPIVOT& _t )
    {
        phred = _t.phred;       // read quality
        ratio = _t.ratio;       // percent identity

        gap = _t.gap;           // number of gaps
        tid = _t.tid;
        odd = _t.odd;           // number of mismatches
        score = _t.score;       // alignment score
        length = _t.length;     // alignment length

        site = _t.site;         // invoke copy constructor

        return( *this );
    }   // operator overloading

    const stPIVOT& operator+=( const stPIVOT& _t )
    {
        phred += _t.phred;      // read quality
        ratio += _t.ratio;      // percent identity

        gap += _t.gap;          // number of gaps
        odd += _t.odd;          // number of mismatches
        score += _t.score;      // alignment score
        length += _t.length;    // alignment length

        for ( unsigned int i = 0; i < _t.site.size(); ++i )
        {
            site.push_back( _t.site[ i ] );
        }   // accumulate the binding sites

        return( *this );
    }   // operator overloading

    double phred;           // read quality
    double ratio;           // percent identity
    unsigned int tid;       // ncbi tid
    unsigned int gap;       // number of gaps
    unsigned int odd;       // number of mismatches
    unsigned int score;     // alignment score
    unsigned int length;    // alignment length
    std::vector<unsigned int> site;
};  // smart container implementation

#endif  // _PIVOT_H
