# Taxonomic Assignment for Whole Metagenome Shotgun (WMGS) Sequencing

Written by [Conrad Shyu](mailto:conradshyu@hotmail.com)<br>
*last revised on July 4, 2020*

> **Note**: This is only the source code repository. The database and sample files are not included because of
their sizes. Instructions are provided to construct the database and translation table, and to download example
WMGS sample files. The repository should include all the necessary files to construct the package.

## Requirements
- PHP (5.2+)
- bowtie2 (the latest)
- GNU C++ compiler (4.4.7+)
- Boost C++ library (1.49+)
- NCBI BlastN (optional; the latest)
- GNU Plot (optional; 4.6+)

## Synopsis
This repository only contains the source code, which implements the weight Shannon equitability index (WSEI) for
the taxonomic assignment of the short reads from whole metagenome shotgun sequencing using the Illumina solid
platform. This approach may be ported to sequencing data generated from other platforms as long as the assumptions
are maintained.

| Filename | Description |
| --- | --- |
| `Makefile` | makefile for the source code |
| `assign.cpp` | taxonomic assignment driver program |
| `samfile.cpp` | bowtie SAM file parser |
| `samfile.h` | header file for bowtie SAM file parser |
| `species.cpp` | taxonomic assignment on the species level |
| `species.h` | header file for the taxnomic assignment program |
| `strain.cpp` | implementation of WSEI |
| `strain.h` | header file for the implementation of WSEI |
| `fas2xlt.cs` | fasta database and translation |
| `README.md` | this file |

To compile the source code, simply type:

`make`

The make program will invoke the Makefile and compiles the source code. Alternatively, to compile the files
manually, issue the command:

```
g++ -I. -O3 samfile.cpp -o samfile -fopenmp
g++ -I. -O3 assign.cpp strain.cpp species.cpp -o assign -fopenmp
```

> Note: The current implementation incorporates automatic multithreading. In other words, the program will
automatically spawn multiple instances based on the number of processor cores. Memory usage is generally not
outrageous. However, WMGS samples are generally quite large and alignment files can be massive.

Assuming that the alignment has been done and output is saved in the SAM format, to perform the analysis, it is
necessary first to parse the alignment SAM file. To parse the SAM file, run the following command:

```
samfile sample.sam sample.summary.csv
```

The parser should requires minimum memory but can be very I/O intensive because it reads and writes files
simultaneously. After the parser completes, run the taxonomic assignment:

```
assign translate.csv sample.summary.csv
```

The taxonomic assignment program will generate two output files, `sample.pivot.csv` and `sample.assign.csv`. The
first file, `sample.pivot.csv`, consolidates the taxonomic assignments on the species level, and the second,
`sample.assign.csv`, lists all candidate taxa that have been identified by the alignment program. Quantitative
statistical analysis should use the first file only. The second file is used to calculate the summary statistics,
i.e., WSEI.

## Construction of Database
Download the database from NCBI ftp server:

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.fna.tar.gz
```
Sequence files are broken into two pieces because of their large sizes. The `bowtie2` databases (complete and
draft genomes) were constructed using the following commands:

- for complete genomes:
```
bowtie2-build -p -f bacteria00.fna bacteria00
bowtie2-build -p -f bacteria01.fna bacteria01
```

- for draft genomes:
```
bowtie2-build -p -f draft00.fna draft00
bowtie2-build -p -f draft01.fna draft01
```

> **Note**: The use of draft genomes is generally not recommended because of potential errors in the sequences. In
addition, the annotation of draft genomes may not be entirely correct either.

> **Note**: It is recommended to remove plasmid sequences before converting to the `bowtie` format.

> **Note**: It is very important to keep both the taxonomy and sequence databases coherent. Otherwise the
assignment program will fail to identify all strains that belong to the same species.

## Map the WMGS Reads
WSEI does not impose any restrictions on how the mapping of short reads is done. The default parameters embedded
in bowtie work just fine. However, it is critically important to map the short reads to the full length genome
database, instead of 16S. The very idea of WSEI is to alleviate the analysis of microbial community from the 16S
doctrine. The use of full length genome also permits more rigorous quantitative statistical analysis on the species
composition and population dynamics.

The default parameters to run bowtie using local alignment:

```
bowtie2 --local --sensitive --threads 4 -x -1 <file1>.fq -2 <file2>.fq -S <file>.sam
```

The parameters can be adjusted to suit specific needs. However, it is not necessary to output all possible hits.
WSEI is only invoked if multiple best hits are present. Otherwise, only the hits with the highest percent identity
is actually used.

Alternatively, it is possible to invoke NCBI BLASTN for the mapping of short reads. BLAST and its variants have
long been the de facto alignment tools. However, they are generally very slow and not suitable for voluminous data
such as WMGS. At this point, the tool does not support alignment output from BLAST yet. This issue will be
addressed on the second release of the tool. The support for BLAST is actually quite simple. It is only necessary
to implement the parser that extract the location of alignment and identity of mapped taxa.

## Example WMGS Files
HMIWGS/HMASM: [Illumina WGS Reads and Assemblies](http://www.hmpdacc.org/HMASM/)

The HMP performed whole metagenomic shotgun sequencing on over 1,200 samples collected from 15 to 18 body sites on
human subjects. The human DNA should have been removed from these samples. WSEI was developed specifically to
analyze the sample using the whole metagenomic shotgun sequencing protocol.

## Author's Comments
The initial implementation of the algorithm requires a slew of other supporting software and files, which can be
very cumbersome to use. The next release will address this issue and, hopefully, simplifies the use of the tool.
Please report any problems or send comments to [me](mailto:conradshyu@hotmail.com).

---
Copyright (C) 2015 [Conrad Shyu](mailto:conradshyu@hotmail.com)<br>
Center for the Study of Biological Complexity<br>
Department of Microbiology and Immunology<br>
Virginia Commonwealth University<br>
Richmond, Virginia 23298<br>

---
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.
