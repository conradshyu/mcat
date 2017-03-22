/*
 * species.cpp
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
 * revised on April 17, 2013
*/

#include <species.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <boost/algorithm/string.hpp>

Species::Species(
    const std::map<unsigned int, stTABLE>& _t,
    const std::map<unsigned int, double>& _w )
{
    const double nMIN = 0.15;
    std::map<unsigned int, stTABLE> t = _t;
    std::map<unsigned int, double> w = _w;
    std::string sid;
    mTaxon.clear(); mIndex.clear();

    for ( std::map<unsigned int, stTABLE>::iterator i = t.begin(); !( i == t.end() ); ++i )
    {
        mTaxon[ ( ( *i ).second ).tid ] = ( ( *i ).second ).species;
    }   // iterate through the records

    for ( std::map<unsigned int, double>::iterator j = w.begin(); !( j == w.end() ); ++j )
    {
        sid = mTaxon.find( ( *j ).first )->second;

        if ( ( *j ).second < nMIN )
        {
            continue;
        }   // eliminate histogram with low index

        mIndex[ sid ] = ( mIndex.find( sid ) == mIndex.end() ) ?
            ( *j ).second : std::max( mIndex.find( sid )->second, ( *j ).second );
    }   // record the weighted shannon index on the species level

    t.clear(); w.clear();
}   // end of copy constructor

Species::~Species()
{
    mTaxon.clear(); mIndex.clear();
}   // default destructor; environmentally conscientious

bool Species::Run( const std::string& _f )
{
    const char* szDELIMIT = ".\n";
    std::string file;
    std::vector<std::string> field;
    mAssign.clear();

    boost::algorithm::split(                // splite the entire string
        field, _f, boost::algorithm::is_any_of( szDELIMIT ) );

    Assign( _f );
    file = field[ 0 ] + ".pivot.csv"; Output( file );
    file = field[ 0 ] + ".assign.csv"; Profile( file );
    mAssign.clear();

    return( true );
}   // end of Run()

/*
 * determine which potential assignment is better
 * histograms have been trimmed to remove taxons that are not actually
 * present in the community
 *
 * 1. assign to taxon if percent identity is higher
 * 2. use shannon index to resolve conflict if ncecessary
 *
 * revised on April 23, 2013
*/
bool Species::Assign(
    const std::string& _r,      // read identification
    const stPIVOT& _p )         // potential assignment
{
    if ( mAssign.find( _r ) == mAssign.end() )
    {
        mAssign[ _r ] = _p; return( true );
    }   // new record has arrived

    unsigned int d1 = mAssign[ _r ].length - mAssign[ _r ].odd;
    unsigned int d2 = _p.length - _p.odd;

    if ( d1 > d2 )
    {
        return( false );
    }   // new assignment is better

    std::string s1 = mTaxon.find( mAssign[ _r ].tid )->second;
    std::string s2 = mTaxon.find( _p.tid )->second;

    if ( mIndex[ s1 ] > mIndex[ s2 ] )
    {
        return( false );
    }   // original shannon index is higher

    mAssign[ _r ] = _p; return( true );
}   // end of Assign()

/*
 * species level assignment
 * summarize the alignment file and generate the output
*/
bool Species::Assign( const std::string& _f )
{
    const double min = 85.0;
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
        std::string rid;
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
            set.tid = static_cast<unsigned int>( ::atoi( field[ 10 ].c_str() ) );   // ncbi tid

            if ( mIndex.find( mTaxon[ set.tid ] ) == mIndex.end() )
            {
                continue;
            }   // histogram not aviable for assignment

            rid = field[ 0 ]; ( set.site ).clear();
            set.ratio = static_cast<double>( ::atof( field[ 1 ].c_str() ) );        // percent identity

            if ( set.ratio < min )
            {
                continue;
            }   // only process good alignment

            set.length = static_cast<unsigned int>( ::atoi( field[ 2 ].c_str() ) ); // alignment length
            set.odd = static_cast<unsigned int>( ::atoi( field[ 3 ].c_str() ) );    // mismatches
            set.gap = static_cast<unsigned int>( ::atoi( field[ 4 ].c_str() ) );    // gaps
            set.phred = static_cast<double>( ::atof( field[ 5 ].c_str() ) );        // read quality
            set.score = static_cast<unsigned int>( ::atoi( field[ 6 ].c_str() ) );  // map quality
            ( set.site ).push_back( set.tid );      // accumulate the count

            #pragma omp critical
            {
                Assign( rid, set );
            }   // the critical region
        } while ( run );    // merge the alignments
    }   // end of the parallel section

    ifs.close(); return( true );
}   // end of Assign()

/*
 * export the assignment of individual read
*/
bool Species::Profile( const std::string& _f)
{
    FILE* of = ::fopen( _f.c_str(), "w" );
    std::string taxon;

    ::fprintf( of, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",       // header
        "Read ID", "Identity", "Alignment Length", "Mismatch", "Gap",
        "Read Quality", "Alignment Quality", "TID", "WSEI", "Taxon" );

    for ( std::map<std::string, stPIVOT>::iterator i = mAssign.begin(); !( i == mAssign.end() ); ++i )
    {
        taxon = mTaxon.find( ( ( *i ).second ).tid )->second;

        ::fprintf( of, "%s,%.2f,%d,%d,%d,%.2f,%d,%d,%.2f,%s\n",
            ( ( *i ).first ).c_str(),   // read identification
            ( ( *i ).second ).ratio,    // average percent identity
            ( ( *i ).second ).length,   // average alignment length
            ( ( *i ).second ).odd,      // average number of mismatches
            ( ( *i ).second ).gap,      // average number of gaps
            ( ( *i ).second ).phred,    // average read quality
            ( ( *i ).second ).score,    // average alignment quality
            ( ( *i ).second ).tid,      // ncbi taxon identification
            mIndex[ taxon ],            // weighted shannon index
            taxon.c_str() );            // species name
    }   // export the assignments

    return( static_cast<bool>( ::fclose( of ) ) );
}   // end of Output()

/*
 * export the contents; just-in-time implementation
*/
bool Species::Output( const std::string& _f )
{
    FILE* of = ::fopen( _f.c_str(), "w" );

    std::map<std::string, stPIVOT> pivot; pivot.clear();
    std::string sid; double count;

    for ( std::map<std::string, stPIVOT>::iterator i = mAssign.begin(); !( i == mAssign.end() ); ++i )
    {
        sid = mTaxon.find( ( ( *i ).second ).tid )->second;
        ( pivot.find( sid ) == pivot.end() ) ?
            pivot[ sid ] = ( *i ).second : pivot[ sid ] += ( *i ).second;
    }   // summarize the assignments first

    ::fprintf( of, "%s,%s,%s,%s,%s,%s,%s,%s\n",     // header
        "Taxon", "Abundance", "Identity", "Alignment Length",
        "Mismatch", "Gap", "Read Quality", "Alignment Quality" );

    for ( std::map<std::string, stPIVOT>::iterator j = pivot.begin(); !( j == pivot.end() ); ++j )
    {
        count = static_cast<double>( ( ( ( *j ).second ).site ).size() );

        ::fprintf( of, "%s,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
            ( ( *j ).first ).c_str(),           // taxon
            static_cast<unsigned int>( ( ( ( *j ).second ).site ).size() ), // abundance
            ( ( *j ).second ).ratio / count,    // average percent identity
            ( ( *j ).second ).length / count,   // average alignment length
            ( ( *j ).second ).odd / count,      // average number of mismatches
            ( ( *j ).second ).gap / count,      // average number of gaps
            ( ( *j ).second ).phred / count,    // average read quality
            ( ( *j ).second ).score / count );  // average alignment quality
    }   // calcualte the weighted shannon index and export the contents

    return( static_cast<bool>( ::fclose( of ) ) );
}   // end of Output()
