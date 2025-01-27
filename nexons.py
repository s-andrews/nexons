#!/usr/bin/env python
import argparse
import tempfile
import subprocess
import os
import sys
import pysam
from progressbar import ProgressBar, Percentage, Bar
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


    #This creates an output directory when it doesn't exist
    outfile_path=Path(options.outfile)
    outfile_path.parent.mkdir(parents=True, exist_ok=True)
    
    
    if options.both_out:
        if options.outfile == "nexons_output.txt":
            gtf_outfile = "nexons_output.gtf"
            custom_outfile = "nexons_output.txt"
        else:
            stripped = options.outfile.replace(".txt", "")
            stripped = stripped.replace(".gtf", "")
            gtf_outfile = stripped + ".gtf"
            custom_outfile = stripped + ".txt"            
        write_gtf_output(quantitations, genes_transcripts_exons, gtf_outfile, options.mincount, splice_info)
        write_output(quantitations, genes_transcripts_exons, custom_outfile, options.mincount, splice_info)

    elif options.gtf_out:
        if options.outfile == "nexons_output.txt":
            outfile = "nexons_output.gtf"
        else:
            outfile = options.outfile
        write_gtf_output(quantitations, genes_transcripts_exons, outfile)
    else:
        write_output(quantitations, genes_transcripts_exons, options.outfile)

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




def write_output(data, gene_annotations, file, mincount, splice_info):
    # The structure for the data is 
    # data[bam_file_name][gene_id][splicing_structure] = {"count":0, "start":[1,2,3], "end":[4,5,6]}
    # 
    # We will have all genes in all BAM files, but might
    # not have all splice forms in all BAM files

    log(f"Writing output to {file} with min count {mincount}")

    bam_files = list(data.keys())
 
    with open(file,"w") as outfile:
        # Write the header
        header = ["File","Gene ID", "Gene Name","Chr","Strand","SplicePattern", "Transcript id","Count","Starts","Ends"]
        outfile.write("\t".join(header))
        outfile.write("\n")

        genes = data[bam_files[0]].keys()

        lines_written = 0

        for gene in genes:
            # We'll make a list of the splices across all samples, and will track the 
            # highest observed value so we can kick out splices which are very 
            # infrequently observed in the whole set.

            splices = {} 

            for bam in bam_files:
                these_splices = data[bam][gene].keys()

                for splice in these_splices:
                    if not splice in splices:
                        splices[splice] = 0

                    if data[bam][gene][splice]["count"] > splices[splice]:
                        splices[splice] = data[bam][gene][splice]["count"]

            # Now extract the subset which pass the filter as we'll report on 
            # all of these

            passed_splices = []

            for splice in splices.keys():
                if splices[splice] >= mincount or options.report_all:
                    passed_splices.append(splice)

            # Now we can go through the splices for all BAM files
            for splice in passed_splices:
                line_values = [
                    gene,
                    gene_annotations[gene]["name"],
                    gene_annotations[gene]["chrom"],
                    gene_annotations[gene]["strand"],
                    ":".join("-".join(str(coord) for coord in junction) for junction in splice),
                    splice_info[gene][splice]["transcript_id"]
                ]

                for bam in bam_files:
                    if splice in data[bam][gene]:
                        splice_line = [
                            str(data[bam][gene][splice]["count"]),
                            ",".join([str(x) for x in data[bam][gene][splice]["start"]]),
                            ",".join([str(x) for x in data[bam][gene][splice]["end"]])
                        ]
                        print("\t".join([bam]+line_values+splice_line), file=outfile)
                
    debug(f"Wrote {len(passed_splices)} splices to {file}")


def write_gtf_output(data, gene_annotations, file, mincount, splice_info):
    # The structure for the data is 
    # data[bam_file_name][gene_id][splicing_structure] = {count:0, start:[1,2,3], end:[4,5,6]}
    # 
    # We will have all genes in all BAM files, but might
    # not have all splice forms in all BAM files

    log(f"Writing GTF output to {file} with min count {mincount}")

    # This copy will be the sum of counts across all samples.
    splice_info_copy = splice_info.copy()
    # set all counts to 0
    for gene in splice_info_copy:
        splice_copy_splices = splice_info_copy[gene].keys()
        for splice_pat in splice_copy_splices: 
            splice_info_copy[gene][splice_pat]["merged_count"] = 0       

    bam_files = list(data.keys())
 
    with open("match_info.txt","w") as report_all_outfile:
        if options.report_all:
            report_all_header = ["seqname", "source", "feature", "start", "end", "id", "exact_count", "merged_count", "splice_pattern"]
            report_all_outfile.write("\t".join(report_all_header))
            report_all_outfile.write("\n")
 
        with open(file,"w") as outfile:
            # Write the header
            header = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
            #header.extend(bam_files)
            outfile.write("\t".join(header))
            outfile.write("\n")

            genes = data[bam_files[0]].keys()

            lines_written = 0

            for gene in genes:
                splices = set() 

                for bam in bam_files:
                    these_splices = data[bam][gene].keys()

                    for splice in these_splices:
                        splices.add(splice)

                gtf_gene_text = "gene_id " + gene

                # Now we can go through the splices for all BAM files
                for splice in splices:
                    #line_values = [gene,gene_annotations[gene]["name"],gene_annotations[gene]["chrom"],gene_annotations[gene]["strand"],splice]
                    splice_start = splice[0][0]
                    splice_end = splice[-1][0]

                    line_values = [gene_annotations[gene]["chrom"], "nexons", "transcript", splice_start, splice_end, 0, gene_annotations[gene]["strand"], 0]

                    splice_text = "splicePattern " + ":".join("-".join(str(coord) for coord in junction) for junction in splice)

                    line_above_min = False
                    for bam in bam_files:
                        if splice in data[bam][gene]:
                           
                            transcript_text = "transcript_id " + str(splice_info[gene][splice]["transcript_id"])
                            
                            debug(f"splice is  {splice}")
                            debug(f"splice info  {splice_info[gene][splice]}")
                            
                            attribute_field = transcript_text + "; " + gtf_gene_text + "; " + splice_text
                            
                            line_values[5] += data[bam][gene][splice]["count"]  # adding the count
                            # if we've got multiple bam files, we want to add up the counts but not keep adding ie. repeating the attribute info.
                            if(len(line_values) == 8):
                            
                                line_values.append(attribute_field)

                            if data[bam][gene][splice]["count"] >= mincount:
                                line_above_min = True
                                
                            # set the count in the splice copy dict
                            splice_info_copy[gene][splice]["merged_count"] += data[bam][gene][splice]["count"]
      
                    if line_above_min:
                        lines_written += 1
                        outfile.write("\t".join([str(x) for x in line_values]))
                        outfile.write("\n")
                        
            if options.report_all:

                splice_info_splices = splice_info_copy[gene].keys()
                for splice_info_splice in splice_info_splices: 
                    splice_is_start = splice_info_splice[0][0]
                    splice_is_end = splice_info_splice[-1][0]
                    other_info = splice_info_copy[gene][splice_info_splice]
                    out_line = [gene, "nexons", "transcript", splice_is_start, splice_is_end, other_info["transcript_id"], other_info["count"], other_info["merged_count"], splice_info_splice]
                    report_all_outfile.write("\t".join([str(x) for x in out_line]))
                    report_all_outfile.write("\n")

    debug(f"Wrote {lines_written} splices to {file}")


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

        if possible_genes:
            breakpoint()

        # for gene in possible_genes:
        #     transcript = gene_matches(exons,gene,flex,endflex)
        #     if transcript is not None:
        #         # Add this to the transcript count
        #         pass


    return counts


def gene_matches(exons,gene,flex,endflex):
    pass


def get_possible_genes(index, chr, start, end):
    # We need to find all genes which completely surround the reported position
    genes = set()
    start_bin = int(start/RESOLUTION)
    end_bin = int(end/RESOLUTION)

    breakpoint()

    for b in range(start_bin,end_bin+1):
        if b >= len(index[chr]):
            # We're past the last gene so don't look any more
            break

        for gene in index[chr][b]:
            if gene["start"]>=start and gene["end"] <= end:
                genes.add(gene)

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
            start = read.reference_start
            current_pos = read.reference_start
            if tuple_operation == 4 or tuple_operation == 5: # soft or hard clip
                start += tuple_length
                current_pos += tuple_length
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

