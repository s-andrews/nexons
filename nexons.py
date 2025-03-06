#!/usr/bin/env python3
import argparse
import sys
import pysam
from pathlib import Path
import json

VERSION = "0.2.devel"

# This is the resolution of the feature indexing - we split the genome
# into bins of this size and use these to quickly find the features we
# need for a given read.
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

    results = []

    # Process each of the BAM files and add their data to the quantitations
    for count,bam_file in enumerate(options.bam):
        log(f"Quantitating {bam_file} ({count+1} of {len(options.bam)})")
        read_lengths,outcomes,quantitations,endflex_observations,innerflex_observations = process_bam_file(genes_transcripts_exons, gene_index, bam_file, options.direction, options.flex, options.endflex)

        results.append(quantitations)

        write_stats_file(bam_file,outcomes, read_lengths,endflex_observations, innerflex_observations, options.outbase)
        write_qc_report(bam_file,outcomes, read_lengths,endflex_observations, innerflex_observations, options.outbase)

        log(f"Summary for {bam_file}:")
        for metric in outcomes:
            log(f"{metric}: {outcomes[metric]}")

    write_output(genes_transcripts_exons,results,options.bam,options.outbase)

def write_stats_file(bam_file, outcomes, read_lengths, endflex, innerflex, outbase):
    outfile = outbase+"_"+bam_file[:-4]+"_stats.txt"
    output = {
        "file": bam_file,
        "outcomes":outcomes,
        "read_lengths": read_lengths,
        "end_flex": endflex,
        "inner_flex": innerflex
    }

    with open(outfile,"wt",encoding="utf8") as out:
        out.write(json.dumps(output))



def write_qc_report(bam_file, outcomes, read_lengths, endflex, innerflex, outbase):
    outfile = outbase+"_"+bam_file[:-4]+"_qc.html"
    template = Path(__file__).parent / "templates/nexons_qc_template.html"

    template_text = ""
    with open(template, "rt", encoding="utf8") as infh:
        for line in infh:
            template_text += line

    template_text = template_text.replace("%%BAMFILE%%",bam_file)

    for metric in outcomes:
        value = f"{outcomes[metric]:,d}"
        template_text = template_text.replace(f"%%{metric}%%",value)
        template_text = template_text.replace(f"%%{metric}_Raw%%",str(outcomes[metric]))


    # We split the read lengths into separate lists for 
    # labels and counts.
    readlengthlabels = []
    readlengthdata = []

    for r in read_lengths:
        readlengthlabels.append(r[0])
        readlengthdata.append(r[1])
    
    template_text = template_text.replace("%%ReadLengthLabels%%",str(readlengthlabels))
    template_text = template_text.replace("%%ReadLengthData%%",str(readlengthdata))


    # We need to assemble the inner and end flexibility data too
    innerflexlabels = []
    innerflexdata = []

    for label,value in innerflex.items():
        innerflexlabels.append(label)
        innerflexdata.append(value)
    
    template_text = template_text.replace("%%InnerFlexLabels%%",str(innerflexlabels))
    template_text = template_text.replace("%%InnerFlexData%%",str(innerflexdata))

    endflexlabels = []
    endflexdata = []

    for label,value in endflex.items():
        endflexlabels.append(label)
        endflexdata.append(value)
    
    template_text = template_text.replace("%%EndFlexLabels%%",str(endflexlabels))
    template_text = template_text.replace("%%EndFlexData%%",str(endflexdata))



    with open(outfile,"wt",encoding="utf8") as out:
        out.write(template_text)


def write_output(genes, quantitations, bam_files, outbase):

    # Quantitation is a list of outputs from the bam files.  Each one of these
    # will be a dictionary with "unique", "partial" and "gene".
    # The unique and partial will have keys of (transcript_id, gene_id) and the
    # gene will just have gene_id.

    for slot in ["unique","partial"]:
        log(f"Processing {slot} output")
        outfile = outbase + "_"+slot+".txt"
        log(f"Outfile is {outfile}")

        # We need a list of all transcript and gene ids used and their detials
        all_ids = set()

        for q in quantitations:
            for id in q[slot]:
                if not id in all_ids:
                    all_ids.add(id)

        # That gets us all of the possible ids.  Now we can go through these in
        # each sample
        
        with open(outfile,"wt",encoding="utf8") as out:
            header = ["Transcript_ID","Gene_ID","Gene_Name","Chr","Start","End","Strand"]
            header.extend(bam_files)

            print("\t".join(header), file=out)

            for id in all_ids:
                gene_id,transcript_id = id
                gene_name = genes[gene_id]["name"]
                chromosome = genes[gene_id]["chrom"]
                start = genes[gene_id]["transcripts"][transcript_id]["start"]
                end = genes[gene_id]["transcripts"][transcript_id]["end"]
                strand = genes[gene_id]["transcripts"][transcript_id]["strand"]

                values = [transcript_id, gene_id,gene_name,chromosome,start,end,strand]

                for q in quantitations:
                    if id in q[slot]:
                        values.append(q[slot][id])
                    else:
                        values.append(0)

                print("\t".join([str(x) for x in values]), file=out)



    for slot in ["gene"]:
        log(f"Processing {slot} output")
        outfile = outbase + "_"+slot+".txt"
        log(f"Outfile is {outfile}")

        # We need a list of gene ids used
        all_ids = set()

        for q in quantitations:
            for id in q[slot]:
                if not id in all_ids:
                    all_ids.add(id)

        # That gets us all of the possible ids.  Now we can go through these in
        # each sample
        
        with open(outfile,"wt",encoding="utf8") as out:
            header = ["Gene_ID","Gene_Name","Chr","Start","End","Strand"]
            header.extend(bam_files)

            print("\t".join(header), file=out)

            for gene_id in all_ids:
                gene_name = genes[gene_id]["name"]
                chromosome = genes[gene_id]["chrom"]
                start = genes[gene_id]["start"]
                end = genes[gene_id]["end"]
                strand = genes[gene_id]["strand"]

                values = [gene_id,gene_name,chromosome,start,end,strand]

                for q in quantitations:
                    if gene_id in q[slot]:
                        values.append(q[slot][gene_id])
                    else:
                        values.append(0)

                print("\t".join([str(x) for x in values]), file=out)


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
    counts = {
        "unique":{},
        "partial":{},
        "gene":{}
    }

    outcomes = {
        "Total_Reads": 0, # Total number of reads in the file
        "No_Alignment" :0, # No reported alignment
        "Primary_Alignment": 0, # Total number of primary alignments
        "Secondary_Alignment": 0, # Total number of secondary alignments
        "No_Gene":0, # Reads not surrounded by a gene
        "No_Hit": 0, # Reads surrounded by gene but not aligning to any transcripts
        "Gene":0, # Reads quantitated at the gene level
        "Multi_Gene":0, # Reads compatible with more than one gene
        "Partial":0, # Reads compatible with one transcript, but not covering all of it
        "Unique":0 # Reads compatible with one transcript and matching completely
    }

    read_lengths = [] # Counts for read lengths for primary alignments in 100bp bins

    # We'll aggregate the observed differences between the exon ends in the reads and
    # the models.
    end_flex_observations = {}
    inner_flex_observations = {}

    for i in range(-flex,flex+1):
        inner_flex_observations[i] = 0

    for i in range(-endflex,endflex+1):
        end_flex_observations[i] = 0


    samfile = pysam.AlignmentFile(bam_file, "rb")

    for read in samfile.fetch(until_eof=True):

        outcomes["Total_Reads"] += 1

        if read.is_unmapped:
            # Nothing to see here
            outcomes["No_Alignment"] += 1
            continue

        if read.is_secondary:
            # This isn't the primary alignment
            outcomes["Secondary_Alignment"] += 1
            continue

        outcomes["Primary_Alignment"] += 1

        length_bin = int(read.query_length / 100)
        if len(read_lengths) < length_bin+1:
            for bin in range(len(read_lengths),length_bin+1):
                read_lengths.append([bin*100,0])

        read_lengths[length_bin][1] += 1

        if not read.reference_name in index:
            # There are no features on this chromsome
            outcomes["No_Gene"] += 1
            continue


        exons = get_exons(read)

        # We need to work out if we want to limit the directionality of the
        # genes we're going to look at.
        # We will allow A (all directions), F (forward only) or R (reverse only) 
        gene_direction="A"

        if direction != "none":
            # We need to care about directionality
            if direction == "same":
                if read.is_reverse:
                    gene_direction = "R"
                else:
                    gene_direction = "F"
            elif direction == "opposite":
                if read.is_reverse:
                    gene_direction = "F"
                else:
                    gene_direction = "R"
            else:
                raise Exception("Unknown direction '"+direction+"'")

        # We want the surrounding genes.  We shorten the read by the amount of 
        # endflex so that we still find genes which surround us if we don't 
        # have perfectly positioned ends.
        possible_genes = get_possible_genes(index, read.reference_name, min(exons[0][0]+endflex, exons[0][1]), max(exons[-1][1]-endflex, exons[-1][0]), gene_direction)

        if not possible_genes:
            outcomes["No_Gene"] += 1
            continue

        found_hit = False
        found_gene_id = None
        found_transcript_id = None
        found_status = None
        best_endflex = None
        best_innerflex = None

        for gene_id in possible_genes:
            # Transcript ID will be the ID of the matched transcript.
            # Status will be one of:
            # 
            # unique, partial, multi 

            transcript_id,status,observed_endflex, observed_innerflex = gene_matches(exons,genes[gene_id],flex,endflex)
            if transcript_id is not None:
                if status=="unique":
                    # We're done here
                    found_hit = True
                    found_gene_id = gene_id
                    found_transcript_id = transcript_id
                    found_status = "unique"
                    best_endflex = observed_endflex
                    best_innerflex = observed_innerflex
                    break

                else:
                    # The logic is the same for both partial and 
                    # multi.  We've locked in the status from this
                    # gene but we could also hit another gene which
                    # would then mean that we don't assign this read
                    # at all.
                    if found_gene_id is not None:
                        # We're done - we can't assign this at all
                        outcomes["Multi_Gene"] += 1
                        found_hit = False
                        break

                    # We can assign this hit but keep looking in case
                    # other genes find anything
                    found_hit = True
                    found_transcript_id = transcript_id
                    found_gene_id = gene_id
                    found_status = status
                    best_endflex = observed_endflex
                    best_innerflex = observed_innerflex
                    continue

         # Now we can increase the appropriate counts
        if not found_hit:
            outcomes["No_Hit"] += 1

        else:
            # There is a hit

            # We can add in the flex values to the total
            for i in best_endflex:
                end_flex_observations[i] += 1

            for i in best_innerflex:
                inner_flex_observations[i] += 1

            # We always increment the gene
            outcomes["Gene"] += 1
            if not found_gene_id in counts["gene"]:
                counts["gene"][found_gene_id] = 1
            else:
                counts["gene"][found_gene_id] += 1

            if found_status == "partial" or found_status == "unique":
                # We increment the partial counts
                outcomes["Partial"] += 1

                if not (found_gene_id,found_transcript_id) in counts["partial"]:
                    counts["partial"][(found_gene_id,found_transcript_id)] = 1
                else:
                    counts["partial"][(found_gene_id,found_transcript_id)] += 1

            if found_status == "unique":
                # We increase the unique count
                outcomes["Unique"] += 1

                if not (found_gene_id,found_transcript_id) in counts["unique"]:
                    counts["unique"][(found_gene_id,found_transcript_id)] = 1
                else:
                    counts["unique"][(found_gene_id,found_transcript_id)] += 1


    # Our read lengths often have high outliers which skew the overall distribution.
    # I'll go through the data until we've hit 5 empty slots - everything above that
    # will be lumped together.
    filtered_read_lengths = []
    zero_count = 0
    stop_adding = False
    for bin in read_lengths:
        if bin[1] == 0:
            zero_count += 1
        else:
            zero_count = 0

        if zero_count == 5:
            stop_adding = True

        if stop_adding:
            filtered_read_lengths[-1][1] += bin[1]
        else:
            filtered_read_lengths.append(bin)


    return (filtered_read_lengths,outcomes,counts,end_flex_observations, inner_flex_observations)


def gene_matches(exons,gene,flex,endflex):

    # We'll look for matches between this structure and 
    # the transcripts in this gene.  We are going to 
    # allow either complete or partial matches. At the
    # end we will return two things - a transcript id
    # and a status
    # 
    # If nothing matches then the transcript will be
    # None
    # 
    # If the read matches the transcript perfectly then
    # the transcript will be the transcript ID and the
    # status will be "unique"
    # 
    # If the read matches partially but only to one 
    # transcript (it doesn't match the others) then
    # the status will be "partial"
    # 
    # If the read partially matches multiple transcripts
    # then the status will be "multi"

    matched_transcript = None
    status = None
    best_endflex = None
    best_innerflex = None

    for transcript in gene["transcripts"].values():

        success, partial, endflex_observed, innerflex_observed = match_exons(exons, transcript["exons"], flex, endflex)

        if success and not partial:
            # It's a full match - we can stop looking
            matched_transcript = transcript["id"]
            status = "unique"
            best_endflex = endflex_observed
            best_innerflex = innerflex_observed
            break

        if success and partial:
            if matched_transcript is None:
                matched_transcript = transcript["id"]
                status = "partial"
                best_endflex = endflex_observed
                best_innerflex = innerflex_observed

            else:
                status = "multi"

    return (matched_transcript,status,best_endflex,best_innerflex)


def match_exons(exons,transcript,flex,endflex):

        # We need to match the exons of this transcript 
        # to the read.

        # We will only allow an exact match, or a sub-match
        # where the read is contained entirely within the 
        # transcript

        # The process will be 
        # 
        # 1. Start from the beginning of the transcript - try to
        # match to the first exon with endflex.
        # 2. If that doesn't work look for an internal match to any 
        # exon using flex.
        # 3. Once we have a match continue it through the exons to
        # see if we can run it to the end.  We'll either get to the
        # end of the transcript (with endflex) or we'll run out within
        # an exon, in which case we'll be a partial match

        full_match = True
        matches = True

        current_transcript_exon = 0
        current_read_exon = 0
        last_transcript_exon = len(transcript)-1
        last_read_exon = len(exons)-1
        end_flex_values = []
        inner_flex_values = []

        while True:

            # We'll step through the different exons of the read to
            # try to match them.

            this_start = exons[current_read_exon][0]
            this_end = exons[current_read_exon][1]

            this_start_flex = flex
            this_end_flex = flex


            # Now we'll add the flexbility which is given to this exon
            if current_transcript_exon == 0:
                this_start_flex = endflex
            
            if current_transcript_exon == last_transcript_exon:
                this_end_flex = endflex

            # See if we match this exon
            start_mismatch =  this_start - transcript[current_transcript_exon][0]
            end_mismatch =  this_end - transcript[current_transcript_exon][1]
            start_matches = abs(start_mismatch) <= this_start_flex
            end_matches = abs(end_mismatch) <= this_end_flex

            if start_matches and end_matches:
                # This exon matches

                # Add the mismatch values to the overall pool
                if current_transcript_exon == 0:
                    end_flex_values.append(start_mismatch)
                else:
                    inner_flex_values.append(start_mismatch)

                            
                if current_transcript_exon == last_transcript_exon:
                    end_flex_values.append(end_mismatch)
                else:
                    inner_flex_values.append(end_mismatch)

                if current_read_exon == last_read_exon:
                    # We're at the end of the match
                    if current_transcript_exon != last_transcript_exon:
                        # We've matched but not to the end of the transcript
                        full_match = False
                    
                    # We don't need to look any further
                    break

                if current_read_exon == 0:
                    # We've matched but are there enough transcript exons
                    # left to accommodate us
                    if len(exons) > (len(transcript)-current_transcript_exon):
                        # This can't possibly work
                        matches = False
                        break

                # We matched and there's more read exons left. 
                current_read_exon += 1
                current_transcript_exon += 1
                continue

            # Not everything matches.

            # We would be in the first exon and not have made the first match yet.
            if current_read_exon == 0 and exons[0][0] > transcript[current_transcript_exon][1]:
                # We can just try the next transcript exon if there is one
                if current_transcript_exon < last_transcript_exon:
                    current_transcript_exon += 1
                    continue
                else:
                    # This isn't going to match
                    matches = False
                    break

            # We could have a partial match.  If we're the first exon of a multi-exon read 
            # then the start could be within the exon.
            if current_read_exon == 0 and len(exons) > 1:
                if end_matches and exons[0][0] >= transcript[current_transcript_exon][0] and exons[0][0] <= transcript[current_transcript_exon][1]:
                    # This could work, but only if we have enough exons left
                    if len(exons) > (len(transcript)-current_transcript_exon):
                        # This can't possibly work
                        matches = False
                        break

                    full_match = False
                    current_read_exon += 1
                    current_transcript_exon += 1
                    continue
                else:
                    matches = False
                    break
                
            # If we're the last exon of a multi-exon read then the end could be within the 
            # exon
            if current_read_exon == last_read_exon and len(exons) > 1:
                if start_matches and exons[current_read_exon][1] >= transcript[current_transcript_exon][0] and exons[current_read_exon][1] <= transcript[current_transcript_exon][1]:
                    full_match = False
                    # This is the last exon so we can stop looking
                    break
                else:
                    matches = False
                    break

                
            # If we're a single exon read then both start and end could be within the exon
            if len(exons) == 1 and exons[0][0] < transcript[current_transcript_exon][1] and exons[0][1] > transcript[current_transcript_exon][0]:
                full_match = False
                break 

            # If we get here then there's no saving us
            matches = False
            break

        return (matches, not full_match, end_flex_values, inner_flex_values)


def get_possible_genes(index, chr, start, end, direction):
    # We need to find all genes which completely surround the reported position
    # Direction is going to be F (forward), R, (reverse) or A (all)
    genes = set()
    start_bin = int(start/RESOLUTION)
    end_bin = int(end/RESOLUTION)

    for b in range(start_bin,end_bin+1):
        if b >= len(index[chr]):
            # We're past the last gene so don't look any more
            break

        for gene in index[chr][b]:
            if direction != "A":
                if gene["strand"] == "+" and direction == "R":
                    continue
                if gene["strand"] == "-" and direction == "F":
                    continue
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

            exon = [start, end]
            

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


            genes[gene_id]["transcripts"][transcript_id]["exons"].append(exon)


    # Before returning the results we need to put the exons
    # for each transcript into order.  We keep everything
    # low - high since that's how the BAM files also structure
    # their results
            
    for gene in genes.values():
        for transcript in gene["transcripts"].values():
            transcript["exons"].sort(key = lambda x: x[0])

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
        "--outbase","-o",
        help="The basename for the output count tables. All outputs will start with this prefix",
        default="./nexons_output"
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
        default=1000, 
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

