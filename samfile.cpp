/*
 * samfile.cpp
 *
 * Written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * Center for the Study of Biological Complexity (CSBC)
 * Department of Microbiology and Immunology
 * Medical College of Virginia
 * Virginia Commonwealth University
 * Richmond, VA 23298
 *
 * simplified api for sam format
 *
 * the store container retain the records:
 * read identification or query tempate name
 * alignment length
 * number of gaps; including deletions and insertions
 * number of mismatches
 * alignment flag
 * 1-base left most mapping position
 * ncbi taxonomy identification
 * ncbi genome identification
 * pred-scale based read quality score
 * mapping quality; alignment quality score
 *
 * to compile:
 * g++ -I. -O3 table.cpp samfile.cpp -o samfile -fopenmp
 * or
 * icc -I. -O2 table.cpp samfile.cpp -o samfile -fopenmp
 *
 * Revised on January 22, 2013
 * Revised on January 24, 2013
 * Revised on February 3, 2013
 * Revised on February 11, 2013
 * Revised on February 13, 2013
*/

#define _DBG_SAMTOOL

#include <omp.h>
#include <table.h>
#include <samfile.h>

#include <map>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

/*
 * default constructor
*/
SamFile::SamFile()
{
}   // default constructor

/*
 * default constructor
*/
SamFile::SamFile(
    const std::map<unsigned int, stTABLE>& _t,  // translation table
    const std::string& _ifs,    // name of alignment file
    const std::string& _ofs )   // name of summary file
{
    Run( _t, _ifs, _ofs );      // multi-threaded version
}   // default constructor

SamFile::~SamFile()
{
}   // default destructor; environmentally conscientious

/*
 * parse the string and assign the variables
 * using gnu regular expression library
 * single thread version
*/
bool SamFile::Run(
    const std::map<unsigned int, stTABLE>& _t,  // translation table
    const std::string& _ifs,            // name of alignment file
    const std::string& _ofs ) const     // name of summary file
{
    const char* szDELIMIT = "\t\n";
    const char* szNCBIBAR = "|";
    const unsigned int nMaxBUFFER = 4096;

    std::ifstream ifs( _ifs.c_str(), std::ios::in );
    FILE* ofs = ::fopen( _ofs.c_str(), "w" );

    ::fprintf( ofs, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
        "Read ID", "Identity", "Length", "Mismatch", "Gaps",
        "Read Quality", "Map Quality", "Left", "Right", "GID", "TID" );

    #pragma omp parallel
    {
        std::vector<std::string> field, ncbi;
        std::map<char, unsigned int> cigar;
        char buffer[ nMaxBUFFER ];
        unsigned int site;
        bool run = false; stSAM sam;

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
            sam.flag = static_cast<unsigned int>( ::atoi( ( field[ 1 ] ).c_str() ) );

            if ( ( field[ 5 ] ).length() < 3 )
            {
                continue;
            }   // not mached properly, according to the aligner

            boost::algorithm::split(                // split the ncbi annotation
                ncbi, field[ 2 ], boost::algorithm::is_any_of( szNCBIBAR ) );
            sam.gid = static_cast<unsigned int>( ::atoi( ( ncbi[ 1 ] ).c_str() ) );

            if ( _t.find( sam.gid ) == _t.end() )
            {
                continue;
            }   // for whatever the reason, gid is not in the table

            ExCIGAR( field[ 5 ], cigar );           // extract cigar string
            sam.qname = field[ 0 ];                 // query template name
            sam.alen = cigar[ 'M' ];                // alignment length
            sam.phred = Sanger( field[ 10 ] );      // phred-scaled score
            sam.off = ExMD( field );                // number of mismatches
            sam.gap = cigar[ 'I' ] + cigar[ 'D' ];  // gaps in alignment
            sam.tid = ( _t.find( sam.gid )->second ).tid;
            site = static_cast<unsigned int>( ::atoi( ( field[ 3 ] ).c_str() ) );
            sam.site = SetBin( site ) + ( _t.find( sam.gid )->second ).start;
            sam.mapq = static_cast<unsigned int>( ::atoi( ( field[ 4 ] ).c_str() ) );
            sam.ratio = static_cast<double>( sam.alen - sam.off + cigar[ 'I' ] )
                / ( sam.alen + cigar[ 'I' ] + cigar[ 'S' ] );

            #pragma omp critical
            {
                ::fprintf( ofs, "\"%s\",%.2f,%d,%d,%d,%.2f,%d,%d,%d,%d,%d\n",
                    sam.qname.c_str(),  // query template name
                    100.0 * sam.ratio,  // percent identity
                    sam.alen,           // alignment length
                    sam.off,            // number of mismatches
                    sam.gap,            // number of gaps; deletions + insertions
                    sam.phred,          // phred-scaled based quality score
                    sam.mapq,           // mapping (alignment) quality
                    sam.site,           // 1-base leftmost mapping position
                    sam.site,           // 1-base rightmost mapping position
                    sam.gid,            // ncbi genome identification
                    sam.tid );          // ncbi taxonomy identification
            }   // the critical region
        } while ( run );    // merge the alignments
    }   // end of the parallel section

    ifs.close(); return( static_cast<bool>( fclose( ofs ) ) );
}   // end of Assign()

/*
 * The SAM FLAGS field, the second field in a SAM record, has multiple bits that
 * describe the paired-end nature of the read and alignment. The first (least
 * significant) bit (1 in decimal, 0x01 in hexidecimal) is set if the read is part
 * of a pair. The second bit (2 in decimal, 0x02 in hexidecimal) is set if the read
 * is part of a pair that aligned in a paired-end fashion. The fourth bit (8 in
 * decimal, 0x08 in hexidecimal) is set if the read is part of a pair and the other
 * mate in the pair had at least one valid alignment. The sixth bit (32 in decimal,
 * 0x20 in hexidecimal) is set if the read is part of a pair and the other mate in
 * the pair aligned to the Crick strand (or, equivalently, if the reverse complement
 * of the other mate aligned to the Watson strand). The seventh bit (64 in decimal,
 * 0x40 in hexidecimal) is set if the read is mate 1 in a pair. The eighth bit (128
 * in decimal, 0x80 in hexidecimal) is set if the read is mate 2 in a pair.
*/
bool SamFile::IsAligned(
    const unsigned int _f ) const
{
    return( ( _f & 0x02 ) ? true : false );
}   // end of IsAlign()

bool SamFile::IsMapped(
    const unsigned int _f ) const
{
    return( ( _f & 0x4 ) ? false : true );
}   // end of IsMapped()

/*
 * if 0x40 and 0x80 are both set, the segment is part of a linear template, but it
 * is neither the first nor the last segment. if both 0x40 and 0x80 are unset, the
 * index of the segment in the template is unknown. this may happen for a non-linear
 * template or the index is lost in data processing.
*/
bool SamFile::IsFirst(
    const unsigned int _f ) const
{
    return( ( _f & 0x40 ) ? true : false );
}   // end of IsFirst()

bool SamFile::IsLast(
    const unsigned int _f ) const
{
    return( ( _f & 0x80 ) ? true : false );
}   // end of IsLast()

/*
 * set the position on the histogram
*/
unsigned int SamFile::SetBin(
    const unsigned int _s ) const
{
    return( static_cast<unsigned int>( ::lrint( _s * 0.001 ) ) );
}   // end of SetBin()

/*
 * calculate the phred-scaled base quality score
*/
double SamFile::Sanger(
    const std::string& _s ) const
{
    const unsigned int nMaxOFFSET = 33;
    double s = 0.0;

    for ( unsigned int i = 0; i < _s.size(); ++i )
    {
        s += ( static_cast<unsigned int>( _s[ i ] ) - nMaxOFFSET );
    }   // accumulate the score

    return( ( _s.size() > 1 ) ? ( s / _s.size() ) : 0.0 );
}   // end of Sanger()

/*
 * extract the cigar string
*/
bool SamFile::ExCIGAR(
    const std::string& _s,
    std::map<char, unsigned int>& _m ) const
{
    char tag;
    const char* cigar = "MIDNSHP=X";
    std::string buffer;
    buffer.clear(); _m.clear();

    for ( unsigned int i = 0; i < ::strlen( cigar ); ++i )
    {
        _m.insert( std::pair<char, unsigned int>( cigar[ i ], 0 ) );
    }   // initialize the map

    for ( unsigned int t = 0; t < _s.size(); ++t )
    {
        tag = _s[ t ];

        if ( ::isalpha( tag ) )
        {
            _m[ tag ] += static_cast<unsigned int>( ::atoi( buffer.c_str() ) );
            buffer.clear(); continue;
        }   // assign the value

        buffer.push_back( tag );    // add a character to the string
    }   // parse the cigar string

    return( true );
}   // end of ExCIGAR()

/*
 * extract the number of mismatches
 *
 * the MD field aims to achieve SNP/indel calling without looking at the
 * reference. for example, a string 10A5^AC6 means from the leftmost
 * reference base in the alignment, there are 10 matches followed by an A on
 * the reference which is different from the aligned read base; the next 5
 * reference bases are matches followed by a 2bp deletion from the reference;
 * the deleted sequence is AC; the last 6 bases are matches. the MD field
 * ought to match the CIGAR string.
*/
unsigned int SamFile::ExMD(
    const std::vector<std::string>& _s ) const
{
    unsigned int m = 0;
    const char* tag = "MD:Z:";

    for ( unsigned int i = 11; i < _s.size(); ++i )
    {
        if ( ( _s[ i ] ).find( tag ) == std::string::npos )
        {
            continue;
        }   // the md tag is not found

        for ( unsigned int t = ::strlen( tag ); t < ( _s[ i ] ).size(); ++t )
        {
            m += ::isalpha( ( _s[ i ] )[ t ] ) ? 1 : 0;
        }   // accumulate the number of mismatches
    }   // optional tages begin at 11

    return( m );
}   // end of ExMD()

#ifdef _DBG_SAMTOOL

/*
 * driver program
 *
 * required parameters:
 * translation table
 * alignment file generated by bowtie
 * ouput filename
*/
int main( int argc, char** argv )
{
    if ( argc < 4 )
    {
        return( 0 );
    }   // check the number of parameters

    const unsigned int nMaxBUFFER = 2048;

    stTABLE a;
    char buffer[ nMaxBUFFER ];
    std::map<unsigned int, stTABLE> table;
    std::ifstream ifs( argv[ 1 ], std::ios::in );

    if ( ifs.fail() )
    {
        return( 0 );
    }   // check the state of stream

    ifs.getline( buffer, nMaxBUFFER );  // skip the header

    while ( ifs.getline( buffer, nMaxBUFFER ) )
    {
        a = buffer; table[ a.gid ] = a;
    }   // parse the file

    ifs.close();
    SamFile s; s.Run( table, argv[ 2 ], argv[ 3 ] );

    return( 0 );
}   // end of main()

#endif  // _DBG_SAMTOOL
