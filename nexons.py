#!/usr/bin/env python
import argparse
import warnings

verbose = False

def main():
    options = get_options()

    global verbose
    verbose = options.verbose

    genes = read_gtf(options.gtf)
    if verbose:
        print(f"Found {len(genes.keys())} genes")

    chromosomes = read_fasta(options.fasta)
    if verbose:
        print(f"Found {len(chromosomes.keys())} chromosomes")


def read_fasta(fasta_file):

    if verbose:
        print(f"Reading sequence from {fasta_file}")

    chromosomes = {}

    with open(fasta_file) as file:

        seqname = None
        sequence = ""

        for line in file:
            if line.startswith(">"):
                if seqname is not None:
                    if seqname in chromosomes:
                        raise Exception(f"Duplicate sequence name {seqname} found in {fasta_file}")

                    chromosomes[seqname] = sequence
                    if verbose:
                        print(f"Added {seqname} {len(sequence)} bp")

                seqname = line.split(" ")[0][1:]
                sequence = ""
            
            else:
                sequence = sequence + line.strip()

        if seqname in chromosomes:
            raise Exception(f"Duplicate sequence name {seqname} found in {fasta_file}")

        chromosomes[seqname] = sequence
        if verbose:
            print(f"Added {seqname} {len(sequence)} bp")

    return chromosomes
        
        



def read_gtf(gtf_file):

    if verbose:
        print(f"Reading genes from {gtf_file}")

    genes = {}

    with open(gtf_file) as file:
        for line in file:
            # Skip comments
            if line.startswith("#"):
                continue

            sections = line.split("\t")

            if len(sections) < 7:
                warnings.warn("Not enough data from line "+line)
                continue

            if sections[2] != "gene":
                continue

            # we can pull out the main information easily enough
            chrom = sections[0]
            start = int(sections[3])
            end = int(sections[4])
            strand = sections[6]

            # For the gene name and gene id we need to delve into the
            # extended comments
            comments = sections[8].split(";")
            gene_id=None
            gene_name=None
            for comment in comments:
                if comment.strip().startswith("gene_id"):
                    gene_id=comment[8:].replace('"','').strip()
                
                if comment.strip().startswith("gene_name"):
                    gene_name=comment.strip()[10:].replace('"','').strip()

            if gene_id is None and gene_name is None:
                warnings.warn(f"No name or id found for gene at {chrom}:{start}-{end}")
                continue

            if gene_id is None:
                warnings.warn(f"Using gene name {gene_name} as ID")
                gene_id = gene_name

            if gene_name is None:
                warnings.warn(f"Using gene ID {gene_id} as name")
                gene_name = gene_id


            genes[gene_id] = {
                "name":gene_name,
                "id": gene_id,
                "chrom": chrom,
                "start": start,
                "end" : end,
                "strand" : strand
            }

    return genes

def get_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "gtf",
        help="A GTF file containing the genes you want to analyse"
    )

    parser.add_argument(
        "fasta",
        help="A multi-fasta file containing the genome sequence"
    )

    parser.add_argument(
        "bam",
        help="One or more BAM files to quantitate",
        nargs="+"
    )

    parser.add_argument(
        "--gene","-g",
        help="The name or ID of a single gene to quantitate"
    )

    parser.add_argument(
        "--verbose","-v",
        help="Produce more verbose output",
        action="store_true"
    )

    return(parser.parse_args())



if __name__ == "__main__":
    main()

