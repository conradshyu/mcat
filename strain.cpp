/*
 * strain.cpp
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
 * revised on April 17, 2013
*/

#include <strain.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

/*
 * constructor
 * set the number of histogram bins
 * set the strain name associated with each tid
*/
Strain::Strain( const std::map<unsigned int, stTABLE>& _t )
{
    std::map<unsigned int, stTABLE> t = _t;
    unsigned int block, tid;

    mBlock.clear(); mTaxon.clear();

    for ( std::map<unsigned int, stTABLE>::iterator i = t.begin(); !( i == t.end() ); ++i )
    {
        tid = ( ( *i ).second ).tid;
        mTaxon[ tid ] = ( ( *i ).second ).strain;
        block = ( ( *i ).second ).end - ( ( *i ).second ).start + 1;

        ( mBlock.find( tid ) == mBlock.end() ) ?
            mBlock[ tid ] = block : mBlock[ tid ] += block; // accumulate genome size
    }   // calcualte the block sizes

    t.clear();      // free up the memory
}   // end of copy constructor

Strain::~Strain()
{
    mBlock.clear(); mTaxon.clear();
}   // default destructor; environmentally conscientious

bool Strain::Run( const std::string& _f )
{
    const char* szDELIMIT = ".\n";
    std::vector<std::string> field;
    mAssign.clear();

    boost::algorithm::split(                // splite the entire string
        field, _f, boost::algorithm::is_any_of( szDELIMIT ) );
    field[ 0 ] += ".strain.csv";

    Assign( _f ); Output( field[ 0 ] );
    mAssign.clear();

    return( true );
}   // end of Run()

/*
 * strain level assignment
 * summarize the alignment file and generate the output
*/
bool Strain::Assign( const std::string& _f )
{
    const char* szDELIMIT = ",\t\n";
    const unsigned int nMaxBUFFER = 2048;
    char header[ nMaxBUFFER ];

    std::ifstream ifs( _f.c_str(), std::ios::in );

    if ( ifs.fail() )
    {
        return( false );
    }   // check the state of stream

    ifs.getline( header, nMaxBUFFER );  // skip the header

    #pragma omp parallel
    {
        char buffer[ nMaxBUFFER ];
        std::vector<std::string> field;
        unsigned int tid;
        stPIVOT set; bool run = false;

        do
        {
            #pragma omp critical
            {
                run = ifs.getline( buffer, nMaxBUFFER );
            }   // the critical region

            if ( !run )
            {
                continue;
            }   // no more data to process

            boost::algorithm::split(                // splite the entire string
                field, buffer, boost::algorithm::is_any_of( szDELIMIT ) );

            ( set.site ).clear();
            set.ratio = static_cast<double>( ::atof( field[ 1 ].c_str() ) );        // percent identity
            set.length = static_cast<unsigned int>( ::atoi( field[ 2 ].c_str() ) ); // alignment length
            set.odd = static_cast<unsigned int>( ::atoi( field[ 3 ].c_str() ) );    // mismatches
            set.gap = static_cast<unsigned int>( ::atoi( field[ 4 ].c_str() ) );    // gaps
            set.phred = static_cast<double>( ::atof( field[ 5 ].c_str() ) );        // read quality
            set.score = static_cast<unsigned int>( ::atoi( field[ 6 ].c_str() ) );  // map quality
            ( set.site ).push_back( static_cast<unsigned int>( ::atoi( field[ 7 ].c_str() ) ) );
            tid = static_cast<unsigned int>( ::atoi( field[ 10 ].c_str() ) );       // ncbi tid

            #pragma omp critical
            {
                ( mAssign.find( tid ) == mAssign.end() ) ?
                    mAssign[ tid ] = set : mAssign[ tid ] += set;
            }   // the critical region
        } while ( run );    // merge the alignments
    }   // end of the parallel section

    ifs.close(); return( true );
}   // end of Assign()

/*
 * export the contents; just-in-time implementation
*/
bool Strain::Output( const std::string& _f )
{
    const double nMIN = 70.0;
    FILE* of = ::fopen( _f.c_str(), "w" );

    mIndex.clear();
    unsigned int tid;
    double weight, optimal, count, block, wsei;

    ::fprintf( of, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",     // header
        "Taxon", "Abundance", "Shannon", "Coverage", "WSEI", "Total Bin", "Identity",
        "Alignment Length", "Mismatch", "Gap", "Read Quality", "Alignment Quality" );

    for ( std::map<unsigned int, stPIVOT>::iterator i = mAssign.begin(); !( i == mAssign.end() ); ++i )
    {
        tid = ( *i ).first;
        block = static_cast<double>( mBlock.find( tid )->second );
        count = static_cast<double>( ( ( ( *i ).second ).site ).size() );
        weight = Weight( ( ( *i ).second ).site, block );
        optimal = ::log( block );
        wsei = Shannon( ( ( *i ).second ).site, weight ) / optimal;

        ::fprintf( of, "%s,%d,%.2f,%.2f,%.2f,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
            ( mTaxon.find( tid )->second ).c_str(),             // taxon
            static_cast<unsigned int>( ( ( ( *i ).second ).site ).size() ), // abundance
            Shannon( ( ( *i ).second ).site, 1.0 ) / optimal,   // conventional shannon
            weight, wsei,                       // coverage and wsei
            mBlock.find( tid )->second,         // total number of bins
            ( ( *i ).second ).ratio / count,    // average percent identity
            ( ( *i ).second ).length / count,   // average alignment length
            ( ( *i ).second ).odd / count,      // average number of mismatches
            ( ( *i ).second ).gap / count,      // average number of gaps
            ( ( *i ).second ).phred / count,    // average read quality
            ( ( *i ).second ).score / count );  // average alignment quality

        if ( ( ( *i ).second ).ratio < ( nMIN * count ) )
        {
            continue;
        }   // only keep the index if percent identity is greate than 85.0

        mIndex[ tid ] = wsei;
    }   // calcualte the weighted shannon index and export the contents

    return( static_cast<bool>( ::fclose( of ) ) );
}   // end of Output()

/*
 * calculate the coverage for a given genome
*/
double Strain::Weight(
    const std::vector<unsigned int>& _s,
    const double _c ) const
{
    std::map<unsigned int, unsigned int> site;

    for ( unsigned int i = 0; i < _s.size(); ++i )
    {
        ( site.find( _s[ i ] ) == site.end() ) ?
            site[ _s[ i ] ] = 1 : site[ _s[ i ] ] += 1;     // accumulate the count
    }   // accumulate the hits for each bin

    return( site.size() / _c );
}   // end of Weight()

/*
 * calculate the conventional/weighted shannon index
 * default weight is 1.0, which is essentially the conventional shannon index
*/
double Strain::Shannon(
    const std::vector<unsigned int>& _s,
    const double _w ) const     // weight; default 1.0
{
    double p, ws = 0.0;
    double t = static_cast<double>( _s.size() );
    std::map<unsigned int, unsigned int> site;

    for ( unsigned int i = 0; i < _s.size(); ++i )
    {
        ( site.find( _s[ i ] ) == site.end() ) ?
            site[ _s[ i ] ] = 1 : site[ _s[ i ] ] += 1;     // accumulate the count
    }   // accumulate the hits for each bin

    for ( std::map<unsigned int, unsigned int>::iterator k = site.begin(); !( k == site.end() ); ++k )
    {
        p = ( *k ).second / t; ws += p * ::log( p );
    }   // caculate the conventional/weighted shannon index

    return( ::fabs( _w * ( ::log( _w ) + ws ) ) );
}   // end of Shannon()

/*
 * retrieve the weighted shannon index
*/
const std::map<unsigned int, double>& Strain::GetIndex() const
{
    return( mIndex );
}   // end of GetIndex()
