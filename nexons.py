#!/usr/bin/env python
import argparse
import warnings
import tempfile
import subprocess
import os

verbose = False

def main():
    options = get_options()

    global verbose
    verbose = options.verbose

    genes = read_gtf(options.gtf, options.gene)
    if verbose:
        print(f"Found {len(genes.keys())} genes")

    chromosomes = read_fasta(options.fasta)
    if verbose:
        print(f"Found {len(chromosomes.keys())} chromosomes")

    quantitations = {}
    for bam_file in options.bam:
        quantitations[bam_file] = process_bam_file(genes, chromosomes, bam_file)


def process_bam_file(genes, chromosomes, bam_file):
    counts = {}

    # Let's keep track of how many reads we've processed
    progress_counter = 0

    # The genes have already been filtered if they're
    # going to be so we can iterate through everything
    for gene_id in genes.keys():

        # For each gene we're going to re-align around
        # the gene region so we'll extract the relevant
        # part of the genome.  We'll start by just using
        # the gene region, but we might want to add some
        # context sequence
        gene = genes[gene_id]

        # Check that we have the sequence for this gene
        if not gene["chrom"] in chromosomes:
            warnings.warn(f"Skipping {gene['name']} as we don't have sequence for chromosome {gene['chrom']}")
            continue

        fasta_file = tempfile.mkstemp(suffix=".fa", dir=".")
        sequence = chromosomes[gene["chrom"]][(gene["start"]-1):gene["end"]]

        with open(fasta_file[1],"w") as out:
            out.write(f">{gene['name']}\n{sequence}\n")

 
        gene_counts = {}
        reads = get_reads(gene,bam_file)

        for read_id in reads.keys():

            # See if we need to print out a progress message
            progress_counter += 1
            if verbose and progress_counter % 10 == 0:
                print("Processed "+str(progress_counter)+" reads")

            segment_string = get_chexons_segment_string(reads[read_id],fasta_file[1],gene["start"])

            if not segment_string in gene_counts:
                gene_counts[segment_string] = 0
            
            gene_counts[segment_string] += 1

        # Clean up the gene sequence file
        os.unlink(fasta_file[1])

        counts[gene["id"]] = gene_counts

    return counts

def get_chexons_segment_string (sequence, genomic_file, position_offset):
    # We need to write the read into a file
    read_file = tempfile.mkstemp(suffix=".fa", dir=".")

    with open(read_file[1],"w") as out:
        out.write(f">read\n{sequence}\n")

    # Now we run chexons to get the data
    chexons_process = subprocess.run(["chexons",read_file[1],genomic_file,"--basename",read_file[1]], check=True, stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    os.unlink(read_file[1]+".comp")

    with open (read_file[1]+".dat") as dat_file:
        locations = [] 
        for line in dat_file:
            if line.startswith("-"):
                continue
            if line.startswith("Seg"):
                continue

            sections = line.split("|")

            if sections[0].strip() == "":
                continue

            if len(sections) < 4:
                continue

            start_end = sections[3].strip().split(" ")
            start = start_end[0]
            end = start_end[-1]

            locations.append([start,end])

    pieces_of_text = []

    for i in range(len(locations)):
        if i==0:
            pieces_of_text.append(str(locations[i][1]))
        elif i==len(locations)-1:
            pieces_of_text.append(str(locations[i][0]))
        else:
            pieces_of_text.append(f"{locations[i][0]}-{locations[i][1]}")

    # Clean up the chexons output
    os.unlink(read_file[1]+".dat")
    os.unlink(read_file[1])


    return ":".join(pieces_of_text)



def get_reads(gene, bam_file):
    # We get all reads which sit within the area of the 
    # bam file.  Since it's a BAM file we need to use 
    # samtools to extract the reads.
    # 
    # We might want to look at doing all of these extractions
    # in a single pass later on.  Doing it one at a time 
    # might be quite inefficient.

    reads = {}

    # The filtered region needs to be in a bed file
    
    bed_file = tempfile.mkstemp(suffix=".bed", dir=".")

    with open(bed_file[1],"w") as out:
        ## TODO: work out how to handle chr names (chr prefix)
        out.write(f"chr{gene['chrom']}\t{gene['start']}\t{gene['end']}\n")

    print(f"Temp file is {bed_file[1]}")

    # Now we can run samtools to get the data
    samtools_process = subprocess.Popen(["samtools","view",bam_file,"-L",bed_file[1]], stdout=subprocess.PIPE)
    
    for line in iter(lambda: samtools_process.stdout.readline(),""):
        sections = line.decode("utf8").split("\t")
        if (len(sections)<9):
            break
        
        if sections[0] in reads:
            warnings.warn("Duplicate read name detected "+sections[0])
            continue

        reads[sections[0]] = sections[9]
    
    # Clean up the bed file
    os.unlink(bed_file[1])
    samtools_process.wait()

    return reads


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
        
        



def read_gtf(gtf_file, gene_filter):

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

            if gene_filter is not None:
                if not (gene_name == gene_filter or gene_id == gene_filter) :
                    continue

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

