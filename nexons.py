#!/usr/bin/env python
import argparse
import tempfile
import subprocess
import os
import sys
import pysam
from pathlib import Path

VERSION = "0.2.devel"
RESOLUTION = 10000



options = argparse.Namespace(verbose=False, quiet=True, report_all=False)

def main():
    global options
    options = get_options()
    
    # Read the details from the GTF for the requested genes
    log(f"Reading GTF {options.gtf}")

    genes_transcripts_exons = read_gtf(options.gtf, options.maxtsl)

    # We should probably build an index so that we can find genes close to reads
    # quickly and efficiently.  This would just be a data structure with a certain
    # resolution which finds all genes within a defined segment of the genome.

    gene_index = build_index(genes_transcripts_exons)

    # Process each of the BAM files and add their data to the quantitations
    quantitations = {}
    for count,bam_file in enumerate(options.bam):
        log(f"Quantitating {bam_file} ({count+1} of {len(options.bam)})")
        quantitations[bam_file] = process_bam_file(genes_transcripts_exons, gene_index, bam_file, options.direction, options.flex, options.endflex)
 
        observations = 0
        for gene in quantitations[bam_file].keys():
            for countdata in quantitations[bam_file][gene].values():
                observations += countdata["count"]

        log(f"Found {observations} valid splices in {bam_file}")


def build_index(genes):
    chr_lengths = {}

    for gene in genes.values():
        if gene["chrom"] not in chr_lengths:
            chr_lengths[gene["chrom"]] = gene["end"]
        else:
            if gene["end"] > chr_lengths[gene["chrom"]]:
                chr_lengths[gene["chrom"]] = gene["end"]


    index = {}

    for chrom in chr_lengths:
        max_bin = int(chr_lengths[chrom]/RESOLUTION)+1
        index[chrom] = []
        for _ in range(max_bin):
            index[chrom].append([])

    # Now we can go through all of the genes appending them to whichever bins they cross
    for gene in genes.values():
        start_bin = int(gene["start"]/RESOLUTION)
        end_bin = int(gene["end"]/RESOLUTION)

        for bin in range(start_bin,end_bin+1):
            index[gene["chrom"]][bin].append(gene)


    return index


def log (message):
    if not options.quiet:
        print("LOG:",message, file=sys.stderr)

def warn (message):
    if not options.suppress_warnings:
        print("WARN:",message, file=sys.stderr)


def debug (message):
    if options.verbose:
        print("DEBUG:",message, file=sys.stderr)


def process_bam_file(genes, index, bam_file, direction, flex, endflex):
    counts = {}

    samfile = pysam.AlignmentFile(bam_file, "rb")

    for read_count,read in enumerate(samfile.fetch(until_eof=True)):

        if not read.reference_name in index:
            # There are no features on this chromsome
            continue

        if read.is_secondary:
            # This isn't the primary alignment
            continue


        exons = get_exons(read)

        possible_genes = get_possible_genes(index, read.reference_name, exons[0][0]-endflex, exons[-1][-1]+endflex)


        for gene_id in possible_genes:
            transcript_id = gene_matches(exons,genes[gene_id],flex,endflex)
            if transcript_id is not None:
                # Add this to the transcript count
                pass


    return counts


def gene_matches(exons,gene,flex,endflex):


    for transcript in gene["transcripts"].values():

        if len([x for x in gene["transcripts"].values()]) == 1 and len(transcript["exons"])>1 and len(transcript["exons"]) == len(exons):
            breakpoint()


        # Check the number of exons matches, since that's really easy
        if not len(transcript["exons"]) == len(exons):
            # Not the same number of exons
            continue


        # Check the ends first as that's easy and uses a different threshold
        min_start=exons[0][0] - endflex
        max_start=min(exons[0][0]+endflex, exons[0][1])

        if not (transcript["start"]>=min_start  and transcript["start"] <= max_start):
            # This can't be a match
            continue


        max_end=exons[-1][1] + endflex
        min_end=max(exons[-1][0], exons[-1][1]-endflex)

        if not (transcript["end"]>=min_end  and transcript["end"] <= max_end):
            # This can't be a match
            continue

        # So we know the ends are OK. 

        # Do the exon positions match well enough
        exons_match = True
        for i in range(len(exons)):
            if i==0:
                # We just check the end
                if not abs(exons[i][1] - transcript["exons"][i][1]) <= flex:
                    exons_match = False
                    break

            elif i==len(exons)-1:
                # We just check the start
                if not abs(exons[i][0] - transcript["exons"][i][0]) <= flex:
                    exons_match = False
                    break
            
            else: 
                # We check start and end
                if not abs(exons[i][1] - transcript["exons"][i][1]) <= flex:
                    exons_match = False
                    break
                if not abs(exons[i][0] - transcript["exons"][i][0]) <= flex:
                    exons_match = False
                    break

        if exons_match:
            breakpoint()
            return transcript["id"]
        
        return None


def get_possible_genes(index, chr, start, end):
    # We need to find all genes which completely surround the reported position
    genes = set()
    start_bin = int(start/RESOLUTION)
    end_bin = int(end/RESOLUTION)

    for b in range(start_bin,end_bin+1):
        if b >= len(index[chr]):
            # We're past the last gene so don't look any more
            break

        for gene in index[chr][b]:
            if start >= gene["start"] and end <= gene["end"]:
                genes.add(gene["id"])

    return genes
    

def get_exons(read):

    exons = []
    start = None
    current_pos = None

    for tuple_operation, tuple_length in read.cigartuples:
        # 0 = match
        # 1 = ins
        # 2 = del
        # 3 = ref_skip
        # 4 = soft_clip
        # 5 = hard_clip
        # 6 = pad
        # 7 = equal
        # 8 = diff
        # 9 = back

        if start is None:
            start = read.reference_start+1 # Reference start is zero based
            current_pos = start
            if tuple_operation == 4 or tuple_operation == 5: # soft or hard clip
                # These aren't included in the reference start position
                # so we just ignore them
                continue
        
        if tuple_operation == 0:
            current_pos += tuple_length
            continue


        if tuple_operation == 1:
            # Read insertion - reference position doesn't change
            continue

        if tuple_operation == 2:
            # Read deletion - increment reference
            current_pos += tuple_length
            continue

        if tuple_operation == 3:
            # Splice operation.  Make a new exon
            exons.append((start,(current_pos-1)))
            current_pos += tuple_length
            start = current_pos
            continue

        else:
            continue

    exons.append((start,(current_pos-1)))

    return exons





def read_gtf(gtf_file, max_tsl):
    debug(f"Reading GTF {gtf_file} with max_tsl {max_tsl}")

    with open(gtf_file) as file:

        genes = {}

        for line in file:

            # Skip comments
            if line.startswith("#"):
                continue

            sections = line.split("\t")

            if len(sections) < 7:
                warn(f"Not enough data from line {line} in {gtf_file}")
                continue

            if sections[2] != "exon":
                continue

            # we can pull out the main information easily enough
            chrom = sections[0]
            start = int(sections[3])
            end = int(sections[4])
            strand = sections[6]

            # For the gene name, gene id and TSL we need to delve into the
            # extended comments
            comments = sections[8].split(";")
            gene_id=None
            gene_name=None
            transcript_id=None
            transcript_name=None
            transcript_support_level=None

            for comment in comments:
                if comment.strip().startswith("gene_id"):
                    gene_id=comment[8:].replace('"','').strip()
                
                if comment.strip().startswith("gene_name"):
                    gene_name=comment.strip()[10:].replace('"','').strip()
                    
                if comment.strip().startswith("transcript_id"):
                    transcript_id=comment[15:].replace('"','').strip()
                          
                if comment.strip().startswith("transcript_name"):
                    transcript_name=comment.strip()[17:].replace('"','').strip()

                if comment.strip().startswith("transcript_support_level"):
                    temp_tsl=comment.strip()[24:].replace('"','').strip().split()[0].strip()
                    if not temp_tsl.isdigit():
                        if temp_tsl=="NA":
                            continue
                        warn(f"Ignoring non-numeric TSL value {temp_tsl}")
                    else:
                        transcript_support_level = int(temp_tsl)

                    
            if gene_id is None and gene_name is None:
                warn(f"No gene name or id found for exon at {chrom}:{start}-{end}")
                continue

            if transcript_id is None and transcript_name is None:
                warn(f"No transcript name or id found for exon at {chrom}:{start}-{end}")
                continue

            if transcript_support_level is None:
                continue

            if transcript_support_level > max_tsl:
                continue

            if gene_id is None:
                gene_id = gene_name

            if gene_name is None:
                gene_name = gene_id
            

            if transcript_id is None:
                transcript_id = transcript_name

            if transcript_name is None:
                transcript_name = transcript_id

            exons = [start, end]
            

            if gene_id not in genes:
                genes[gene_id] = {
                    "name": gene_name,
                    "id": gene_id,
                    "chrom": chrom,
                    "start": start,
                    "end" : end,
                    "strand": strand,
                    "transcripts": {
                    }
                }
            else:
                if start < genes[gene_id]["start"]:
                    genes[gene_id]["start"] = start
                if end > genes[gene_id]["end"]:
                    genes[gene_id]["end"] = end


            if transcript_id not in genes[gene_id]["transcripts"]:
                genes[gene_id]["transcripts"][transcript_id] = {
                    "name": transcript_name,
                    "id": transcript_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "exons" : []
                }
            else:
                if start < genes[gene_id]["transcripts"][transcript_id]["start"]:
                    genes[gene_id]["transcripts"][transcript_id]["start"] = start

                if end > genes[gene_id]["transcripts"][transcript_id]["end"]:
                    genes[gene_id]["transcripts"][transcript_id]["end"] = end


            genes[gene_id]["transcripts"][transcript_id]["exons"].append([start,end])


    return genes




def get_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "gtf",
        help="A GTF file containing the genes you want to analyse"
    )

    parser.add_argument(
        "--maxtsl",
        help="Maximum transcript support level to analyse",
        type=int,
        default=1
    )

    parser.add_argument(
        "bam",
        help="One or more BAM files to quantitate",
        nargs="+"
    )

    parser.add_argument(
        "--outfile","-o",
        help="The file to write the output count table to",
        default="nexons_output.txt"
    )

    parser.add_argument(
        "--flex","-f",
        help="How many bases different can exon boundaries be and still merge them",
        default=10, 
        type=int
    )

    parser.add_argument(
        "--endflex","-e",
        help="How many bases different can transcript ends be and still merge them",
        default=50, 
        type=int
    )

    parser.add_argument(
        "--mincount","-m",
        help="What is the minimum number of observations for any given variant to report",
        default=2, 
        type=int
    )

    parser.add_argument(
        "--direction","-d",
        help="The directionality of the library (none, same, opposing)",
        default="none"
    )

    parser.add_argument(
        "--verbose","-v",
        help="Produce more verbose output",
        action="store_true"
    )

    parser.add_argument(
        "--quiet",
        help="Suppress all messages",
        action="store_true"
    )

    parser.add_argument(
        "--suppress_warnings",
        help="Suppress warnings (eg about lack of names or ids)",
        action="store_true"
    )

    parser.add_argument(
        "--version",
        help="Print version and exit",
        action="version",
        version=f"Nexons version {VERSION}"
    )

    return(parser.parse_args())



if __name__ == "__main__":
    main()

