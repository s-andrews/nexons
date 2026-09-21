![Nexons Logo](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/nexons_logo_path.svg)

# Introduction
Nexons is a program to quantitate RNA-Seq data from nanopore sequencing runs.  It takes in a BAM file aligned with a suitable spliced aligner, and a GTF file of gene annotations from the genome to which the BAM file was aligned and creates a series of count tables from transcript or gene level matches with differing degress of confidence. It also generates QC reports to summarise the findings from each file.

![Read Fate](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/read_fate.png)
![Alignment Fate](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/alignment_fate.png)
![Alignment Directionality](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/alignment_directionality.png)
![Read Length and Coverage](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/read_length_and_coverage.png)
![Exon Flex](https://raw.githubusercontent.com/s-andrews/nexons/refs/heads/master/images/exon_flex.png)


After running nexons you can use [nexons viewer](https://www.bioinformatics.babraham.ac.uk/nexons-viewer/) to look at the results.


# Installation
Nexons is a python program which should work under any recent version of python3.  It depends on the following non-core packages

* pysam (https://github.com/pysam-developers/pysam)

# Usage
The basic usage of the program is simply:

```
nexons.py [gtf_file] [bam1] [bam2]...
```

This will quantitate the bam files into output files starting with ```nexons_output```

Use `--parallel 4` to process up to four BAM files concurrently (default: 1).
Each file is analysed in its own worker process using the existing analysis.
Per-file reports and aggregated tables are written as usual, with aggregate
columns retaining the input BAM order. Input BAMs must have distinct basenames
when using parallel processing, because output names are based on those names.
Each worker holds annotation data, so higher concurrency uses more memory.

## Additional options

```
usage: nexons.py [-h] [--parallel PARALLEL] [--maxtsl MAXTSL] [--outbase OUTBASE] [--flex FLEX]
                 [--endflex ENDFLEX] [--direction DIRECTION] [--allqc] [--noannotate] [--verbose] [--quiet]
                 [--suppress_warnings] [--version]
                 gtf bam [bam ...]

positional arguments:
  gtf                   A GTF file containing the genes you want to analyse
  bam                   One or more BAM files to quantitate

options:
  -h, --help            show this help message and exit
  --parallel PARALLEL   Number of BAM files to process concurrently (default 1)
  --maxtsl MAXTSL       Maximum transcript support level to analyse (default 2)
  --outbase OUTBASE, -o OUTBASE
                        The basename for the output count tables. All outputs will start with this prefix (default
                        ./nexons_output)
  --flex FLEX, -f FLEX  How many bases different can exon boundaries be and still merge them (default 3)
  --endflex ENDFLEX, -e ENDFLEX
                        How many bases different can transcript ends be and still merge them (default 5000)
  --direction DIRECTION, -d DIRECTION
                        The directionality of the library [none, same, opposing] (default none)
  --allqc               Write individual sample QC HTML reports as well as the combined report
  --noannotate          Skip the production of annotated BAM files
  --verbose, -v         Produce more verbose output
  --quiet               Suppress all messages
  --suppress_warnings   Suppress warnings (eg about lack of names or ids)
  --version             Print version and exit
```

# Output files

By default, the program writes a combined HTML QC report and per-sample text
QC files. Add `--allqc` to also write individual HTML QC reports for every sample.

## `<outbase>_combined_qc.html`
An additional QC report comparing all samples from the run, in input order.
It includes every outcome count, stacked horizontal charts for read fate and directionality, five separate
alignment-class bar charts, and a line per sample for read lengths,
transcript coverage and exon flexibility. Each bar chart has a counts/percentage control; percentages use all BAM entries
for that sample. Alignment classes retain their overlapping counts. Like the individual
reports, charts load Chart.js from the CDN and need an internet connection.

## ```nexons_output_unique.txt```
A count table of unique hits which span the full length of a transcript in the GTF file

## ```nexons_output_partial.txt```
A count table of unique full and partial hits to a transcript in the GTF file.  This includeds all of the counts in the ```nexons_output_unique.txt``` file as well as the unambiguously mapped partial hits.

## ```nexons_output_gene.txt```
A gene level count table including all hits where a read matches part of a transcript, or multiple transcripts, within the same gene.  Hits which map to more than one possible gene are not included.

## ```nexons_output_[filename]_qc.html```
Created only when `--allqc` is supplied.
An HTML QC report summarising the matches found in each file

### GTF annotation statistics

Run `nexons_gtf_stats.py` to produce a three-column TSV (`Metric`, `Value`,
`Count`) from a GTF or gzipped GTF:

```bash
python3 nexons_gtf_stats.py annotation.gtf.gz --maxtsl 2 --outbase annotation_stats
python3 nexons_gtf_stats.py annotation.gtf.gz --chromosome chr1 --outbase chr1_stats
```

The output prefix is set with `--outbase` (or `-o`), defaulting to
`nexons_gtf_stats`. Each run writes `<outbase>.txt` (tab-delimited metrics) and
`<outbase>.html` (an interactive report). The report includes the input file,
command and effective options, totals, biotype doughnut charts, and distribution
line charts. Drag, Ctrl-scroll or pinch to zoom on each line chart’s x axis;
Shift-drag pans, and Reset zoom restores the full range. Like the other Nexons
reports, the HTML loads Bootstrap and Chart.js from a CDN and requires an
internet connection for chart libraries (including the zoom plugin).
TSL filtering matches `nexons.py`: the default maximum is 2, missing or `NA`
TSLs are excluded, and MANE_Select, Ensembl_Canonical, gencode_primary and
 gencode_basic annotations are treated as TSL 1. `--maxtsl 0` disables filtering.
Filtering applies to exon records, as in the main script. Counts include only
transcripts with retained exons and genes containing those transcripts.
`--chromosome` (also `--chrom`) requires the exact chromosome name and stops at
the next chromosome after the selected block; the input must be grouped by
chromosome.

Gene/transcript counts include biotypes from `gene_biotype`/`gene_type` and
`transcript_biotype`/`transcript_type`, including attributes on gene and transcript
records. Missing transcript biotypes fall back to the gene biotype; otherwise
missing biotypes are labelled `unknown`.

Exon lengths include both endpoints. Mature transcript lengths use 100 bp bins
labelled `0-99bp`, `100-199bp`, etc. Pair comparisons are unordered, within a gene
and strand, with first/last exons determined in transcript direction. Internal
alternate splice comparisons exclude first/last exons and identical ends.
Terminal comparisons require multi-exon transcripts and include identical
boundaries in the zero bin. Each qualifying transcript/exon pair contributes;
shared events across multiple pairs are counted multiple times.

Alternate splice lengths use single-base bins, ending in `100+bp`. Transcript
start lengths use `0bp`, `1-5bp`, `6-10bp`, etc., ending in `996-999bp`
and `1000+bp`. Transcript end lengths use `0bp`, `1-50bp`, `51-100bp`, etc.,
ending in `9951-9999bp` and `10000+bp`. Both terminal distributions keep
exact zero differences separate; the HTML Exclude zero toggle hides that bin. Final categories include the
cap itself and all larger distances. These three distance histograms include
zero-count bins; other distributions report observed bins only.

Use `--distinct` to count each unordered pair of distinct
terminal positions once per gene and strand in the start/end distributions.
Only positions represented by multi-exon transcripts with a shared first-exon
end (for starts) or last-exon start (for ends), in transcript direction, are
compared. A position pair qualifying through multiple shared boundaries is
still counted once. Different position pairs with the same distance remain
separate observations. Repeated identical positions collapse to one position,
so this mode has no zero-distance pairs. Other metrics are unchanged; by default,
all qualifying transcript pairs are counted. Both text and HTML outputs use the
selected mode, which is also recorded in the HTML summary.

```bash
python3 nexons_gtf_stats.py annotation.gtf.gz --distinct --outbase distinct_stats
```
