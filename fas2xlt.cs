/*
 * fas2xlt.cs
 *
 * C# implementatoin that splits FASTA files and remove plasmid. The program
 * will prepare the full-length genomes for taxonomic assignment. After
 * download the full length genome databse from NCBI, unzip and run the program
 * to recursively traverse into each directory and process individual sequence.
 * There is no need to manually merge all the sequences first. This program will
 * also generate the translation table required by the algorithm. It is
 * important that all the following files are on the same directory: names.dmp,
 * nodes.dmp, gi_taxid_nucl.dmp, and the program.
 *
 * Download locations:
 * ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
 * ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.zip
 * ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.fna.tar.gz
 *
 * Example usage:
 * mono fsa2xlt.exe bacteria
 *
 * This command will generate a FASTA file named, bacteria.fna, and a
 * translation table, bacteria.csv.
 *
 * Note: This is memory efficient implementation. The memory usage does not
 * depend on the size of database.
 *
 * Note: NCBI taxonomy files are incredibly large. It will take some time for
 * the program to complete.
 *
 * Copyright (C) 2015   Conrad Shyu (shyu4751@yahoo.com)
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
*/

public class xlt
{
    private System.Int32 td;    // taxonomic identification
    private System.Int32 rk;    // taxonomic rank
    private System.UInt32 sz;   // size of the genome
    private System.String st;   // strain name
    private System.String sp;   // species name

    public xlt( System.String _k, System.UInt32 _s )
    {
        sp = ( ( ( ( _k.Split( '|' ) )[ 4 ] ).Split( ',' ) )[ 0 ] ).Trim();
        td = -1; rk = -1; st = sp; sz = _s;
    }   // class constructor

    public System.Int32 tid
    {
        get { return( td ); } set { td = value; }
    }   // set or get the taxonomic identification

    public System.Int32 rank
    {
        get { return( rk ); } set { rk = value; }
    }   // set or get the species rank

    public System.UInt32 length
    {
        get { return( sz ); } set { sz = value; }
    }   // set or get the genome length

    public System.String strain
    {
        get { return( st ); } set { st = value; }
    }   // set or get the strain name

    public System.String species
    {
        get { return( sp ); } set { sp = value; }
    }   // set or get the species name
}   // translation table structure

public class fas2xlt
{
    private const System.String NCBI_NAMES = "names.dmp";
    private const System.String NCBI_NODES = "nodes.dmp";
    private const System.String NCBI_TAXID = "gi_taxid_nucl.dmp";
    private const System.Double WMGS_BLOCK = 1000;

    static System.UInt32 set_tid(
        ref System.Collections.Generic.Dictionary<System.Int32, xlt> _x,
        System.String _f = NCBI_TAXID )
    {
        System.UInt32 i = 0;

        try
        {
            using ( System.IO.StreamReader r = new System.IO.StreamReader( _f ) )
            {
                System.String b;

                for ( i = 0; ( b = r.ReadLine() ) != null; ++i )
                {
                    if ( ( b.Trim() ).Length < 2 )
                    {
                        continue;
                    }   // skip empty line

                    System.Int32 key = System.Convert.ToInt32( ( ( b.Split( '\t' ) )[ 0 ] ).Trim() );

                    if ( !_x.ContainsKey( key ) )
                    {
                        continue;
                    }   // not in the record

                    ( _x[ key ] ).tid = System.Convert.ToInt32( ( ( b.Split( '\t' ) )[ 1 ] ).Trim() );
                    ( _x[ key ] ).rank = ( _x[ key ] ).tid;
                }   // set the taxonomic identification and rank
            }
        }
        catch ( System.Exception e )
        {
            System.Console.WriteLine( "Exception: {0}", e.Message );
        }   // an exception has occurred

        return( i );    // total number of records processed
    }   // set the taxonomic identification

    static System.UInt32 set_rank(
        ref System.Collections.Generic.Dictionary<System.Int32, xlt> _x,
        System.String _f = NCBI_NODES )
    {
        System.UInt32 i = 0;
        System.Collections.Generic.Dictionary<System.Int32, System.Int32> tid =
            new System.Collections.Generic.Dictionary<System.Int32, System.Int32>();

        foreach( System.Collections.Generic.KeyValuePair<System.Int32, xlt> a in _x )
        {
            if ( tid.ContainsKey( ( a.Value ).tid ) )
            {
                continue;
            }   // duplicate taxonomic identification

            tid.Add( ( a.Value ).tid, -1 );
        }   // populate the taxonomic identification

        try
        {
            using ( System.IO.StreamReader r = new System.IO.StreamReader( _f ) )
            {
                System.String b;

                for ( i = 0; ( b = r.ReadLine() ) != null; ++i )
                {
                    if ( ( b.Trim() ).Length < 2 )
                    {
                        continue;
                    }   // skip empty line

                    System.Int32 key = System.Convert.ToInt32( ( ( b.Split( '|' ) )[ 0 ] ).Trim() );

                    if ( !tid.ContainsKey( key ) )
                    {
                        continue;
                    }   // not in the record

                    tid[ key ] = ( ( ( b.Split( '|' ) )[ 2 ] ).Trim() ).Contains( "no rank" ) ?
                        System.Convert.ToInt32( ( ( b.Split( '|' ) )[ 1 ] ).Trim() ) : key;
                }   // set the taxonomic identification and rank
            }
        }
        catch ( System.Exception e )
        {
            System.Console.WriteLine( "Exception: {0}", e.Message );
        }   // an exception has occurred

        foreach( System.Collections.Generic.KeyValuePair<System.Int32, xlt> c in _x )
        {
            ( c.Value ).rank = tid[ ( c.Value ).tid ];
        }   // populate the taxonomic identification

        return( i );    // total number of records processed
    }   // set the taxonomic rank

    static System.UInt32 set_name(
        ref System.Collections.Generic.Dictionary<System.Int32, xlt> _x,
        System.String _f = NCBI_NAMES )
    {
        System.UInt32 i = 0;
        System.Collections.Generic.Dictionary<System.Int32, System.String> rank =
            new System.Collections.Generic.Dictionary<System.Int32, System.String>();

        foreach( System.Collections.Generic.KeyValuePair<System.Int32, xlt> a in _x )
        {
            if ( rank.ContainsKey( ( a.Value ).rank ) )
            {
                continue;
            }   // duplicate taxonomic identification

            rank.Add( ( a.Value ).rank, ( a.Value ).species );
        }   // populate the taxonomic identification

        try
        {
            using ( System.IO.StreamReader r = new System.IO.StreamReader( _f ) )
            {
                System.String b; System.String[] p;

                for ( i = 0; ( b = r.ReadLine() ) != null; ++i )
                {
                    p = b.Split( '|' );

                    if ( p.Length < 4 )
                    {
                        continue;
                    }   // skip empty line

                    if ( !( p[ 3 ] ).Contains( "scientific" ) )
                    {
                        continue;
                    }   // wrong field

                    if ( rank.ContainsKey( System.Convert.ToInt32( p[ 0 ].Trim() ) ) )
                    {
                        rank[ System.Convert.ToInt32( p[ 0 ].Trim() ) ] = p[ 1 ].Trim();
                    }   // assign the name
                }   // set the taxonomic identification and rank
            }
        }
        catch ( System.Exception e )
        {
            System.Console.WriteLine( "Exception: {0}", e.Message );
        }   // an exception has occurred

        foreach( System.Collections.Generic.KeyValuePair<System.Int32, xlt> c in _x )
        {
            ( c.Value ).species = rank[ ( c.Value ).rank ];
        }   // populate the species name

        return( i );    // total number of records processed
    }   // set the species name

    static System.Int32 fasta(
        System.Collections.Generic.Dictionary<System.String, System.String> _l,
        System.String _f = "complete.fna" )
    {
        System.String k = "";
        System.Text.StringBuilder s = new System.Text.StringBuilder();
        System.String[] text = System.IO.File.ReadAllLines( _f );

        foreach ( System.String a in text )
        {
            if ( ( a.Trim() ).Length < 2 )
            {
                continue;
            }   // skip the empty line

            if ( !( a[ 0 ] == '>' ) )
            {
                s.Append( a.Trim() ); continue;      
            }   // accumulate sequence

            k = a.Trim( new System.Char[] { '>', ' ' } );
        }   // process the contents in the array

        _l.Add( k, ( s.ToString() ).Trim() ); return( _l.Count );
    }   // read the fasta file

    public static int Main( System.String[] args )
    {
        if ( args.Length < 1 )
        {
            System.Console.WriteLine( "require: output" ); return( 0 );
        }   // check the required parameters

        System.String of = args[ 0 ] + ".fna";
        System.String ot = args[ 0 ] + ".csv";

        System.IO.DirectoryInfo[] dir = new System.IO.DirectoryInfo( @"." ).GetDirectories();
        System.Collections.Generic.Dictionary<System.Int32, xlt> csv =
            new System.Collections.Generic.Dictionary<System.Int32, xlt>();

        try
        {
            using ( System.IO.StreamWriter w = new System.IO.StreamWriter( of, false ) )
            {
                foreach( System.IO.DirectoryInfo d in dir )
                {
                    System.Collections.Generic.Dictionary<System.String, System.String> list =
                        new System.Collections.Generic.Dictionary<System.String, System.String>();

                    System.Console.WriteLine( "directory: {0}", d.Name );
                    System.String[] fna = System.IO.Directory.GetFiles(
                        d.Name, "*.fna", System.IO.SearchOption.AllDirectories );

                    foreach( System.String f in fna )
                    {
                        System.Console.Write( "file: {0} ", f.Split( new System.Char[] { '\\', '/' } )[ 1 ] );
                        fasta( list, f ); System.Console.WriteLine( "completed" );
                    }   // process each file

                    foreach( System.Collections.Generic.KeyValuePair<System.String, System.String> a in list )
                    {
                        if ( ( a.Key ).Contains( "plasmid" ) )
                        {
                            continue;
                        }   // skip plasmid sequences

                        w.Write( ">{0}\n{1}\n", a.Key, a.Value );
                        csv.Add( System.Convert.ToInt32( ( ( a.Key ).Split( '|' ) )[ 1 ] ),
                            new xlt( a.Key, ( System.UInt32 )( ( a.Value ).Length ) ) );
                    }   // write the fasta file
                }   // process each directory
            }   // export the sequence into a file
        }
        catch ( System.Exception e )
        {
            System.Console.WriteLine( "Exception: {0}", e.Message );
        }   // an exception has occurred

        System.Console.Write( "processing {0} ... ", NCBI_TAXID );
        set_tid( ref csv ); System.Console.WriteLine( "completed" );

        System.Console.Write( "processing {0} ... ", NCBI_NODES );
        set_rank( ref csv ); System.Console.WriteLine( "completed" );

        System.Console.Write( "processing {0} ... ", NCBI_NAMES );
        set_name( ref csv ); System.Console.WriteLine( "completed" );

        try
        {
            using ( System.IO.StreamWriter w = new System.IO.StreamWriter( ot, false ) )
            {
                w.WriteLine( "\"GID\",\"TID\",\"Size\",\"Start\",\"End\",\"Strain\",\"Species\"" );

                foreach( System.Collections.Generic.KeyValuePair<System.Int32, xlt> i in csv )
                {
                    if ( ( i.Value ).tid < 0 )
                    {
                        continue;
                    }   // skip incomplete record

                    w.WriteLine( "{0},{1},{2},0,{3},\"{4}\",\"{5}\"",
                        i.Key, ( i.Value ).tid, ( i.Value ).length,
                        System.Math.Round( ( i.Value ).length / WMGS_BLOCK ),
                        ( i.Value ).strain, ( i.Value ).species );
                }   // dump the contents
            }   // write the translation table
        }
        catch ( System.Exception e )
        {
            System.Console.WriteLine( "Exception: {0}", e.Message );
        }   // an exception has occurred

        return( 1 );
    }   // main procedure
}
