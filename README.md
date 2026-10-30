# Taxonomic Assignment for Whole Metagenome Shotgun Sequencing

Written by [Conrad Shyu](mailto:conradshyu@hotmail.com)

This repository contains C++ tools for processing Bowtie2 SAM alignments and assigning metagenomic reads to bacterial strains and species. The strain-level analysis uses the weighted Shannon equitability index (WSEI) to help resolve ambiguous alignments.

The sequence and taxonomy databases, translation table, and example sequencing data are not included. They must be prepared separately.

## Requirements

- GNU Make
- A C++ compiler with OpenMP support (the Makefile uses `g++` and `-fopenmp`)
- Boost C++ headers (the source includes Boost.Algorithm)
- Bowtie2 to build an index and align reads
- Mono or a .NET compiler/runtime, only if using the optional `fas2xlt.cs` database-preparation helper

BLAST output is not supported; the parser expects Bowtie2-compatible SAM records.

## Build

From the repository directory, run:

```sh
make
```

This builds two programs:

- `samfile` parses a SAM alignment file and writes a summary CSV.
- `assign` calculates strain-level WSEI and species-level assignments from that summary.

To build without Make:

```sh
g++ -I. -O3 samfile.cpp -o samfile -fopenmp
g++ -I. -O3 assign.cpp strain.cpp species.cpp -o assign -fopenmp
```

## Prepare the database

The optional `fas2xlt.cs` helper creates a FASTA file and the translation CSV required by the analysis tools. It expects `names.dmp`, `nodes.dmp`, and `gi_taxid_nucl.dmp` in the working directory, alongside directories containing the extracted genome FASTA files. It searches for `.fna` files under those directories and accepts an output prefix:

```sh
mcs fas2xlt.cs -out:fas2xlt.exe
mono fas2xlt.exe bacteria
```

The command writes `bacteria.fna` and `bacteria.csv`. Build a Bowtie2 index from the FASTA file:

```sh
bowtie2-build --threads 4 -f bacteria.fna bacteria
```

Keep the translation table and sequence index in sync. The table maps the genome identifiers found in the alignment references to NCBI taxonomic identifiers and names.

## Align and analyze reads

For paired-end reads, map against the full-length genome index (not a 16S database). The `--no-hd` option omits SAM header lines, which this parser does not handle:

```sh
bowtie2 --local --sensitive --threads 4 --no-hd \
  -x bacteria -1 reads_1.fq -2 reads_2.fq -S sample.sam
```

Parse the alignments, then run the assignment:

```sh
samfile bacteria.csv sample.sam sample.summary.csv
assign bacteria.csv sample.summary.csv
```

The analysis produces three CSV files using the summary filename's prefix:

- `sample.strain.csv` contains strain-level abundance, coverage, Shannon index, and WSEI.
- `sample.pivot.csv` aggregates read assignments and alignment statistics by species.
- `sample.assign.csv` lists the selected species assignment and alignment statistics for each read.

Use the species-level `*.pivot.csv` file for quantitative species-abundance analysis. The per-read `*.assign.csv` file provides assignment details; `*.strain.csv` contains strain-level WSEI statistics.

## Source files

| File | Purpose |
| --- | --- |
| `Makefile` | Builds the C++ programs |
| `samfile.cpp`, `samfile.h` | Parse SAM alignments and produce a summary CSV |
| `assign.cpp` | Load the translation table and run both assignment stages |
| `strain.cpp`, `strain.h` | Calculate strain-level coverage and WSEI |
| `species.cpp`, `species.h` | Resolve assignments and summarize results by species |
| `table.h` | Translation-table record definition and CSV parsing |
| `pivot.h` | Alignment and assignment summary record |
| `fas2xlt.cs` | Optional FASTA and translation-table preparation helper |
| `LICENSE` | License terms |

## Example data

Example human microbiome whole-genome shotgun reads and assemblies are available from the [Human Microbiome Project](http://www.hmpdacc.org/HMASM/). The project removed human DNA from its samples; WSEI was developed for whole-metagenome shotgun data.

For questions or problems, contact [Conrad Shyu](mailto:conradshyu@hotmail.com).
