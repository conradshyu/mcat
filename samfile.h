/*
 * samfile.h
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
 * Revised on January 22, 2013
 * Revised on January 24, 2013
 * Revised on February 3, 2013
 * Revised on February 14, 2013
*/

#ifndef _SAMFILE_H
#define _SAMFILE_H

#include <table.h>

#include <map>
#include <list>
#include <vector>
#include <string>

struct stSAM
{
    stSAM()
    {
        qname.clear();
    }   // default constructor

    stSAM( const stSAM& _s )
    {
        *this = _s;
    }   // copy constructor

    ~stSAM()
    {
        qname.clear();
    }   // default destructor; environmentally conscientious

    std::string qname;      // rid; query template name
    unsigned int alen;      // alignment length
    unsigned int gap;       // number of gaps; deletions and insertions
    unsigned int off;       // number of mismatches
    unsigned int flag;      // alignment flag reported by bowtie
    unsigned int site;      // 1-base leftmost mapping position
    unsigned int tid;       // ncbi taxonomy identification
    unsigned int gid;       // ncbi genome identification
    unsigned int mapq;      // mapping quality; alignment quality score
    double ratio;           // percent identity
    double phred;           // phred-scale based read quality score
};  // definition of container class

/*
 * interface class for the container
*/
class SamFile
{
public:
    SamFile();
    SamFile(
        const std::map<unsigned int, stTABLE>&,
        const std::string&, const std::string& );
    ~SamFile();

    bool Run(
        const std::map<unsigned int, stTABLE>&,
        const std::string&, const std::string& ) const;

private:
    /*
     * A typical use of a function object is in writing callback functions.
     * A callback in procedural languages, such as C, may be performed by using
     * function pointers. However it can be difficult or awkward to pass a
     * state into or out of the callback function. This restriction also
     * inhibits more dynamic behavior of the function. A function object solves
     * those problems since the function is really a façade for a full object,
     * carrying its own state.
    */
    struct SortEx
    {
        bool TestTID( const stSAM& _a, const stSAM& _b ) const
        {
            return( _a.tid > _b.tid );
        }   // end of TestTID()

        bool TestSize( const stSAM& _a, const stSAM& _b ) const
        {
            unsigned int i = _a.alen - _a.off;
            unsigned int j = _b.alen - _b.off;

            return( ( i > j ) ? true : false );
        }   // end of TestSize()

        bool TestName( const stSAM& _a, const stSAM& _b ) const
        {
            return( ( _a.qname.compare( _b.qname ) == 0 ) ? true : false );
        }   // end of TestName()

        /*
         * sorting criteria
         * 1. read identification
         * 2. ncbi taxonomy identification
         * 3. alignmnet lengths (without mismatches)
        */
        bool operator()( const stSAM& _a, const stSAM& _b ) const
        {
            if ( TestName( _a, _b ) )
            {
                return( true );
            }   // check the query name

            if ( TestTID( _a, _b ) )
            {
                return( true );
            }   // check the taxnomy identification

            return( TestSize( _a, _b ) );
        }   // end of operator overloading
    };  // end of class SortEx

    double Sanger( const std::string& ) const;
    unsigned int SetBin( const unsigned int ) const;
    unsigned int ExMD( const std::vector<std::string>& ) const;

    bool IsLast( const unsigned int ) const;
    bool IsFirst( const unsigned int ) const;
    bool IsMapped( const unsigned int ) const;
    bool IsAligned( const unsigned int ) const;
    bool ExCIGAR( const std::string&, std::map<char, unsigned int>& ) const;
};  // end of class definition

#endif  // _SAMTOOL_H
