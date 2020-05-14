# nexons
A pipeline for quantitating transcript level abundances from nanopore sequence data

This code is a wrapper around the chexons program which maps cDNA sequences to genomic segments.  In this pipeline we use chexons to perform a local re-alignment of nanopore sequencing reads which have been genomically mapped with a general purpose mapper to extract specific splice junction coordinates.  We then collate and quantitate these both within and between samples.

Quick Start
-----------

```
./nexons.py \
--gene=Pkm \
--verbose \
test_data/Mus_musculus.GRCm38.99.gtf \
test_data/Mus_musculus.GRCm38.dna.chromosome.9.fa \
test_data/test.bam
```
