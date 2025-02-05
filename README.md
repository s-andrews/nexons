![Nexons Logo](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/nexons_logo_path.svg)

# Introduction
Nexons is a program to quantitate RNA-Seq data from nanopore sequencing runs.  It takes in a BAM file aligned with a suitable spliced aligner, and a GTF file of gene annotations from the genome to which the BAM file was aligned and creates a series of count tables from transcript or gene level matches with differing degress of confidence. It also generates QC reports to summarise the findings from each file.

![Read Fate](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/read_fate.png)
![Alignment Fate](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/alignment_fate.png)

# Installation
Nexons is a python program which should work under any recent version of python3.  It depends on the following non-core packages

* pysam (https://github.com/pysam-developers/pysam)

# Usage
The basic usage of the program is simply:

```
nexons.py [gtf_file] [bam1] [bam2]...
```

This will quantitate the bam files into output files starting with ```nexons_output```

## Additional options

```
usage: nexons.py [-h] [--maxtsl MAXTSL] [--outbase OUTBASE] [--flex FLEX] [--endflex ENDFLEX] 
[--direction DIRECTION] [--verbose] [--quiet] [--suppress_warnings] [--version] gtf bam [bam ...]

positional arguments:
  gtf                   A GTF file containing the genes you want to analyse
  bam                   One or more BAM files to quantitate

options:
  -h, --help            show this help message and exit

  --maxtsl MAXTSL       Maximum transcript support level to analyse

  --outbase -o          The basename for the output count tables. All outputs will start with this prefix

  --flex -f             How many bases different can exon boundaries be and still merge them

  --endflex -e          How many bases different can transcript ends be and still merge them
  --direction -d        The directionality of the library (none, same, opposing)
  --verbose, -v         Produce more verbose output
  --quiet               Suppress all messages
  --version             Print version and exit
```

# Output files

## ```nexons_output_unique.txt```
A count table of unique hits which span the full length of a transcript in the GTF file

## ```nexons_output_partial.txt```
A count table of unique full and partial hits to a transcript in the GTF file.  This includeds all of the counts in the ```nexons_output_unique.txt``` file as well as the unambiguously mapped partial hits.

## ```nexons_output_gene.txt```
A gene level count table including all hits where a read matches part of a transcript, or multiple transcripts, within the same gene.  Hits which map to more than one possible gene are not included.

## ```nexons_output_[filename]_qc.html```
An HTML QC report summarising the matches found in each file


