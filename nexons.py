#!/usr/bin/env python
import argparse
import tempfile
import subprocess
import os
import sys
from progressbar import ProgressBar, Percentage, Bar
import pickle

VERSION = "0.1.devel"

options = argparse.Namespace(verbose=False, quiet=True, report_all=False)

def main():
    global options
    options = get_options()
    
    # Read the details from the GTF for the requested genes
    # for now, we're going to read the gtf multiple times, as first we want the genes, then transcripts, then exons. 
    # If the gtf was ordered, we could do it in one pass, but safer to do it 3 times - might be slow though...
    log(f"Reading genes from {options.gtf} with gene list {options.gene}")

    genes_transcripts_exons = read_gtf(options.gtf, options.gene)

    # now all the exons have been imported from the gtf, we can convert the splice patterns to a compatible format.
    debug("Converting splice patterns")
    for gene_id in genes_transcripts_exons:
        for transcript_id in genes_transcripts_exons[gene_id]["transcripts"]:
            exon_set = genes_transcripts_exons[gene_id]["transcripts"][transcript_id]["exons"]
            splices = convert_splice_pattern(exon_set)
            genes_transcripts_exons[gene_id]["transcripts"][transcript_id]["splice_patterns"] = splices
     
    log(f"Found {len(genes_transcripts_exons.keys())} genes to quantitate")
    
    if len(genes_transcripts_exons.keys())==0:
        raise Exception("No genes found")

    log(f"Reading chromosomes from {options.fasta}")
    chromosomes = read_fasta(options.fasta)
    
    log(f"Found {len(chromosomes.keys())} chromosomes")

    if len(chromosomes.keys())==0:
        raise Exception("No chromosomes found")

    # Process each of the BAM files and add their data to the quantitations
    quantitations = {}
    for count,bam_file in enumerate(options.bam):
        log(f"Quantitating {bam_file} ({count+1} of {len(options.bam)})")
        quantitations[bam_file] = process_bam_file(genes_transcripts_exons, chromosomes, bam_file, options.direction, options.minexons, options.mincoverage, options.mapthreshold)
 
        observations = 0
        for gene in quantitations[bam_file].keys():
            for countdata in quantitations[bam_file][gene].values():
                observations += countdata["count"]

        log(f"Found {observations} valid splices in {bam_file}")

    with open("quantitations.pkl","wb") as out:
        pickle.dump(quantitations,out)

    with open("genes.pkl","wb") as out:
        pickle.dump(genes_transcripts_exons,out)


    log("Collating splice variants")

    quantiations, splice_info = collate_splice_variants(quantitations,options.flexibility, genes_transcripts_exons)
    
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
        write_gtf_output(quantitations, genes_transcripts_exons, outfile, options.mincount, splice_info)
    else:
        write_output(quantitations, genes_transcripts_exons, options.outfile, options.mincount, splice_info)

def log (message):
    if not options.quiet:
        print("LOG:",message, file=sys.stderr)

def warn (message):
    if not options.suppress_warnings:
        print("WARN:",message, file=sys.stderr)


def debug (message):
    if options.verbose:
        print("DEBUG:",message, file=sys.stderr)



def collate_splice_variants(data, flexibility, genes_transcripts_exons):
    # The structure for the data is 
    # data[bam_file_name][gene_id][splicing_structure_tuple] = {count:0, start:1, end:2}
    # 
    # We will have all genes in all BAM files, but might
    # not have all splice forms in all BAM files
    # 
    # Here we aim to produce a reduced set of splice variants
    # by combining variants which differ by only a few positions

    debug("Merging similar variants")

    bam_files = list(data.keys())

    genes = data[bam_files[0]].keys()

    # Build the structure for the merged data
    merged_data = {}
    for bam in data.keys():
        merged_data[bam] = {}
        for gene in genes:
            merged_data[bam][gene] = {}

    splice_counts = {} 

    # Now work our way through each gene doing the merging
    for gene in genes:
       
        splice_counts[gene] = {}
       
        # simplifying the splice counts - genes_transcripts_exons[gene]["transcripts"] contains a load of info for that transcript - a dict of start, stop, exons, splice_patterns (maybe it doesn't need that but we'll leave it for now).
       
        # seed with the transcripts from the gtf file
        for transcript_id in genes_transcripts_exons[gene]["transcripts"]:

            # structure is splice_counts[gene][splice_pattern]
            splice_counts[gene][genes_transcripts_exons[gene]["transcripts"][transcript_id]["splice_patterns"]] = {
                "transcript_id": transcript_id,
                "count": 0,
                "strand": genes_transcripts_exons[gene]["transcripts"][transcript_id]["strand"]
            }
            
        #print(f"strand = {genes_transcripts_exons[gene]['transcripts'][transcript_id]['strand']}")
        #print("\n ============ splice_counts[gene].keys()==================")
        #print(splice_counts[gene].keys())
        #print(splice_counts[gene].keys())
        #try saving them separately, sorting, then adding them in
       
        # get all the splice_patterns
        for bam in bam_files:
            these_splices = data[bam][gene].keys()
           # print(f"data[bam][gene]: {data[bam][gene]}")
            
            for splice in these_splices:
                
                if not splice in splice_counts[gene].keys():
                    splice_counts[gene][splice] = {
                        "transcript_id": "unknown",#_transcript,
                        "count": 0,
                        "strand": "tbc"
                    }
                    #unknown_transcript += 1               
                splice_counts[gene][splice]["count"] += data[bam][gene][splice]["count"]
                
              #  print(f"splice counts info: {splice_counts[gene][splice]}")
                
        # split the dictionary so that we order the known transcripts separately from the novel ones.
        known_splices = {}
        novel_splices = {}
        
        for splice in splice_counts[gene].keys():
            if splice_counts[gene][splice]["transcript_id"] == "unknown":
                novel_splices[splice] = splice_counts[gene][splice]["count"]
            else:
                known_splices[splice] = splice_counts[gene][splice]["count"]

        #print("\  novel_splice_counts")
        #print(novel_splices)
        
        # these become sorted lists
        known_splices_sorted = sorted(known_splices.keys(), key=lambda x:known_splices[x], reverse=True)
        novel_splices_sorted = sorted(novel_splices.keys(), key=lambda x:novel_splices[x], reverse=True)
        
        # join them back together
        all_splices = known_splices_sorted + novel_splices_sorted
        
        #print("\  novel_splice_counts_sorted")
        #print(novel_splices_sorted)
        
        
        # Now we can collate the splice counts.  We'll get
        # back a hash of the original name and the name we're
        # going to use.
        # 
        # We'll give the function the ordered set of splices so
        # it always uses the most frequently observed one if 
        # there is duplication -
        
        splice_name_map = create_splice_name_map(all_splices, flexibility) 

        # From this map we can now build up a new merged set
        # of quantitations which put together the similar splices
        for bam in bam_files:
            these_splices = data[bam][gene].keys()

            for splice in these_splices:
                used_splice = splice_name_map[splice]
                if not used_splice in merged_data[bam][gene]:
                    merged_data[bam][gene][used_splice] = {"count":0, "start":[], "end":[]}
                merged_data[bam][gene][used_splice]["count"] += data[bam][gene][splice]["count"]
                merged_data[bam][gene][used_splice]["start"].extend(data[bam][gene][splice]["start"])
                merged_data[bam][gene][used_splice]["end"].extend(data[bam][gene][splice]["end"])

    # We should also sort the start and end lists
    for bam in merged_data.keys():
        for gene in merged_data[bam].keys():
            for splice in merged_data[bam][gene].keys():
                merged_data[bam][gene][splice]["start"].sort()
                merged_data[bam][gene][splice]["end"].sort()

    return [merged_data, splice_counts]



def create_splice_name_map(splices, flexibility):
    debug(f"Merging {len(splices)} different splice sites")
    # This takes an ordered list of splice tuples and matches
    # them on the basis of how similar they are.  We work 
    # our way down the list trying to match to previous splices

    # This is what we'll give back.  Every string will be in this
    # and we'll match it either to itself or a more popular 
    # equivalent string


    map_to_return = {}

    # These are the ordered set of patterns which we retained
    # ie we didn't merge them with another pattern.
    valid_splice_patterns = []

    for splice_to_test in splices:

        # Now we try to map the parsed segment to existing segments

        # We now test the existing splice patterns to see if they're
        # a good enough match to this one.
        found_hit = False

        # If there aren't the same number of exons then it's defintitely
        # not a match.
        for test_segment in valid_splice_patterns:
            if len(test_segment) != len(splice_to_test):
                continue

            # Try to match the coordinates of the segments
            too_far = False
            # Iterate though the segments
            for i in range(len(splice_to_test)):
                # Iterate through the start/end positions
                for j in range(len(splice_to_test[i])):
                    if abs(splice_to_test[i][j]-test_segment[i][j]) > flexibility:
                        too_far = True
                        break
                if too_far:
                    break
            if not too_far:
                # We can use this as a match
                #print(f"Merged:\n{splice}\ninto\n{test_segment['string']}\n\n")
                map_to_return[splice_to_test] = test_segment
                found_hit = True
                break
        
        if not found_hit:
            # Nothing was close enough, so enter this as a new reference
            map_to_return[splice_to_test] = splice_to_test
            valid_splice_patterns.append(splice_to_test)

    if options.verbose:    
        debug(f"Produced {len(valid_splice_patterns)} deduplicated splices")
    
    return map_to_return



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
                    str(splice),
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
    # data[bam_file_name][gene_id][splicing_structure] = count
    # 
    # We will have all genes in all BAM files, but might
    # not have all splice forms in all BAM files

    log(f"Writing GTF output to {file} with min count {mincount}")

    # make a copy so we can get get counts, including 0s.
    # The counts in the original are the number of splices that mapped exactly, not taking into account the flexibility factor.
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
                    splice_split = splice.split(":")
                    splice_start = splice_split[0]
                    splice_end = splice_split.pop()

                    line_values = [gene_annotations[gene]["chrom"], "nexons", "transcript", splice_start, splice_end, 0, gene_annotations[gene]["strand"], 0]

                    splice_text = "splicePattern " + splice

                    line_above_min = False
                    for bam in bam_files:
                        if splice in data[bam][gene]:
                           
                            transcript_text = "transcript_id " + str(splice_info[gene][splice]["transcript_id"])
                            
                            debug(f"splice is  {splice}")
                            debug(f"splice info  {splice_info[gene][splice]}")
                            
                            attribute_field = transcript_text + "; " + gtf_gene_text + "; " + splice_text
                            
                            line_values[5] += data[bam][gene][splice]  # adding the count
                            # if we've got multiple bam files, we want to add up the counts but not keep adding ie. repeating the attribute info.
                            if(len(line_values) == 8):
                            
                                line_values.append(attribute_field)

                            if data[bam][gene][splice] >= mincount:
                                line_above_min = True
                                
                            # set the count in the splice copy dict
                            splice_info_copy[gene][splice]["merged_count"] += data[bam][gene][splice]
      
                    if line_above_min:
                        lines_written += 1
                        outfile.write("\t".join([str(x) for x in line_values]))
                        outfile.write("\n")
                        
            if options.report_all:

                splice_info_splices = splice_info_copy[gene].keys()
                for splice_info_splice in splice_info_splices: 
                    splice_is_split = splice_info_splice.split(":")
                    splice_is_start = splice_is_split[0]
                    splice_is_end = splice_is_split.pop()
                    other_info = splice_info_copy[gene][splice_info_splice]
                    out_line = [gene, "nexons", "transcript", splice_is_start, splice_is_end, other_info["transcript_id"], other_info["count"], other_info["merged_count"], splice_info_splice]
                    report_all_outfile.write("\t".join([str(x) for x in out_line]))
                    report_all_outfile.write("\n")

    debug(f"Wrote {lines_written} splices to {file}")


def process_bam_file(genes, chromosomes, bam_file, direction, min_exons, min_coverage, map_threshold):
    counts = {}

    # we need a set of splices that are passed in - that need to be in a dictionary of gene_ids as there may be 
    # multiple genes we're looking for.

    # The genes have already been filtered if they're
    # going to be so we can iterate through everything
    for gene_id in genes.keys():

        # For each gene we're going to re-align around
        # the gene region so we'll extract the relevant
        # part of the genome.  We'll start by just using
        # the gene region, but we might want to add some
        # context sequence
        gene = genes[gene_id]

        log(f"Quantitating {gene_id} in {bam_file}")

        # Check that we have the sequence for this gene
        if not gene["chrom"] in chromosomes:
            warn(f"Skipping {gene['name']} as we don't have sequence for chromosome {gene['chrom']}")
            continue

        fasta_file = tempfile.mkstemp(suffix=".fa", dir="/dev/shm")
        sequence = chromosomes[gene["chrom"]][(gene["start"]-1):gene["end"]]

        with open(fasta_file[1],"w") as out:
            out.write(f">{gene['name']}\n{sequence}\n")
 

        # We're going to keep track of the set of results from chexons.  The 
        # gene counts data structure has a key of a tuple of the splice boundaries
        # and a value which is a dict with a count of number of times observed, as
        # well as lists of the observed start and end positions so we can later
        # match up transcripts.
        gene_counts = {}  
        reads = get_reads(gene,bam_file,direction)

        log(f"Found {len(reads)} reads for gene {gene_id} in {bam_file}")

        # Let's keep track of how many reads we've processed
        progress_counter = 0

        # Let's count the failures
        reasons_for_failure = {}

        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=len(reads))
        if not options.quiet:
            pbar.start()

        for read_id in reads.keys():

            # See if we need to print out a progress message
            progress_counter += 1
            if not options.quiet:
                pbar.update(progress_counter)
            
            # Running chexons will fail if there are no matches between 
            # the read and the gene (or too short to make a match segment)
            # We therefore need to catch this and skip any reads which 
            # trigger it

            try:
                chexons_result = get_chexons_segment_string(reads[read_id],fasta_file[1],gene, min_exons, min_coverage, map_threshold)
                
                
                # This will either be a dict with the match details in if it worked
                # or it will be a string starting with "FAIL" if not.  Yes, I've heard
                # about exceptions, but this works right now.
                if type(chexons_result) is dict:
                    splice_boundaries = chexons_result["splice_boundaries"]
                    if not splice_boundaries in gene_counts:
                        gene_counts[splice_boundaries] = {"count": 0, "start":[], "end": []}
                
                    gene_counts[splice_boundaries]["count"] += 1
                    gene_counts[splice_boundaries]["start"].append(chexons_result["start"])
                    gene_counts[splice_boundaries]["end"].append(chexons_result["end"])

                else:
                    # It's a failure string
                    if not chexons_result in reasons_for_failure:
                        reasons_for_failure[chexons_result] = 0
                    
                    reasons_for_failure[chexons_result] += 1
                            
            except Exception as e:
                warn(f"[WARNING] Chexons failed with {e}")

        if not options.quiet:
            pbar.finish()

        log("Reasons for chexons result rejection")
        for message,count in reasons_for_failure.items():
            log(f"{count} {message}")

        # Clean up the gene sequence file
        os.unlink(fasta_file[1])
        
        counts[gene["id"]] = gene_counts

    return counts

def get_chexons_segment_string (sequence, genomic_file, gene, min_exons, min_coverage, map_threshold):

    position_offset = gene["start"]
    # We need to write the read into a file
    read_file = tempfile.mkstemp(suffix=".fa", dir="/dev/shm")

    with os.fdopen(read_file[0],"w") as out:
        out.write(f">read\n{sequence}\n")


    # Now we run chexons to get the data
    try:
        chexons_process = subprocess.run(["chexons",read_file[1],genomic_file,"--basename",read_file[1]], check=True, stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    except Exception as ex:
        # If it's failed we need to clean up anything left behind
        os.unlink(read_file[1]+".comp")
        os.unlink(read_file[1]+".dat")
        os.remove(read_file[1])
        raise ex


    os.unlink(read_file[1]+".comp")

    count = 0
    start_cDNA = 0
    end_cDNA = 0
    cDNA_length = 0
    full_sequence_length = len(sequence)

    with open (read_file[1]+".dat") as dat_file:
        locations = [] 
        for line in dat_file:
            #print(line, end=None)
        
            if line.startswith("-"):
                continue
            if line.startswith("Seg"):
                count = 0
                continue

            sections = line.split("|")
            
           # print(f"chexons sections: {sections}")
            
            if sections[0].strip() == "":
                continue

            if len(sections) < 4:
                continue
           
            #if sections[2].strip().split(" ")[1] == "F":
            #    potential_direction = "F"
            #else: 
                #potential_direction = sections[2].strip().split(" ")[1]
            #    print(f"potential_direction: {sections[2]}")

            count += 1

            if count == 1:
                start_end_cDNA = sections[1].strip().split(" ")
                start_cDNA = int(start_end_cDNA[0])

            start_end = sections[3].strip().split(" ")
            start = int(start_end[0])+position_offset-1
            end = int(start_end[-1])+position_offset-1

            locations.append([start,end])
        
            start_end_cDNA = sections[1].strip().split(" ")           
            end_cDNA = int(start_end_cDNA[-1])
        cDNA_length = end_cDNA - start_cDNA
    
    # Clean up the chexons output
    os.unlink(read_file[1]+".dat")
    os.remove(read_file[1])

    #print(f"cDNA length = {cDNA_length}, full sequence length = {full_sequence_length}")
    proportion_mapped = cDNA_length/full_sequence_length
    
    if options.verbose_proportions:
        log(f"Proportion of full sequence mapped = {proportion_mapped:.3f}")

    # Check that enough of the sequence has been mapped
    if proportion_mapped < map_threshold:
        if options.verbose_proportions:
            log(f"Discarding sequence as proportion mapped is too low {proportion_mapped} vs {map_threshold}")
        return "FAIL: Transcript coverage too low"

    # Check that we've got enough coverage
    if abs(locations[0][0] - locations[-1][-1]) / abs(gene["end"]-gene["start"]) < min_coverage:
        debug(f"Discarding sequence as coverage of gene was too low {abs(locations[0][0] - locations[-1][-1])} vs {min_coverage}")
        return "FAIL: Gene coverage too low" 

    # Check that we've got enough exons to keep this
    if len(locations) < min_exons:
        debug(f"Discarding sequence as number of exons found {len(locations)} is lower than the threshold {min_exons}")
        return "FAIL: Not enough exons"


    # Our splice boundaries are going to be the locations without the start and end position
    splice_boundaries = []

    for i in range(len(locations)):
        if i==0:
            splice_boundaries.append((locations[i][1],))
        elif i==len(locations)-1:
            splice_boundaries.append((locations[i][0],))
        else:
            splice_boundaries.append((locations[i][0],locations[i][1]))


    return_data = {
        "start":locations[0][0],
        "end": locations[-1][1],
        "splice_boundaries": tuple(splice_boundaries)
    }

    return return_data



def get_reads(gene, bam_file, direction):
    # We get all reads which sit within the area of the 
    # bam file.  Since it's a BAM file we need to use 
    # samtools to extract the reads.
    # 
    # We might want to look at doing all of these extractions
    # in a single pass later on.  Doing it one at a time 
    # might be quite inefficient.

    reads = {}

    # The filtered region needs to be in a bed file
    
    bed_file = tempfile.mkstemp(suffix=".bed", dir="/dev/shm")

    with open(bed_file[1],"w") as out:
    
        if(options.no_chr_prefix):
            out.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\n") 
            
        else:
            # Some chromosome names will already start with chr but we need
            # to add it here if it doesn't already
            if gene['chrom'].startswith("chr"):
                out.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\n")

            else:
                out.write(f"chr{gene['chrom']}\t{gene['start']}\t{gene['end']}\n")


    # Now we can run samtools to get the data.  We require that the
    # read overlaps the relevant region, but we also check the direction
    # of the read.  Different library preps generate different directions
    # so this can be 'none', 'same' or 'opposing'.
    #
    # The filter for forward strand reads is -f 16 and reverse is -F 16
    # the library is opposing strand specific so if we have a forward
    # strand gene we want -f and we'd use -F for reverse strand.

    if direction == "none":
        #print("Launching samtools with","samtools","view",bam_file,"-L",bed_file[1])
        samtools_process = subprocess.Popen(["samtools","view",bam_file,"-L",bed_file[1]], stdout=subprocess.PIPE)

    elif direction == "opposing":
        strand_filter_string = "-f" if gene["strand"] == "+" else "-F"
        debug("Launching samtools with","samtools","view",bam_file,"-L",bed_file[1],strand_filter_string,"16")
        samtools_process = subprocess.Popen(["samtools","view",bam_file,"-L",bed_file[1],strand_filter_string,"16"], stdout=subprocess.PIPE)

    elif direction == "same":
        strand_filter_string = "-F" if gene["strand"] == "+" else "-f"
        debug("Launching samtools with","samtools","view",bam_file,"-L",bed_file[1],strand_filter_string,"16")
        samtools_process = subprocess.Popen(["samtools","view",bam_file,"-L",bed_file[1],strand_filter_string,"16"], stdout=subprocess.PIPE)

    else:
        raise Exception("Didn't understand directionality "+direction)

    
    for line in iter(lambda: samtools_process.stdout.readline(),""):
        sections = line.decode("utf8").split("\t")
        
        if (len(sections)<9):
            debug(f"Only {len(sections)} sections in {line} from samtools, so exiting")
            break
       # print(sections[3])
        if sections[0] in reads:
            warn(f"Duplicate read name {sections[0]} detected")
            continue

        reads[sections[0]] = sections[9]
    
    # Clean up the bed file
    os.unlink(bed_file[1])
    samtools_process.wait()

    return reads


def read_fasta(fasta_file):

    debug(f"Reading sequence from {fasta_file}")

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
                    debug(f"Added {seqname} {len(sequence)} bp")

                seqname = line.split(" ")[0][1:]
                seqname = seqname.rstrip() # Laura - removing new line character in fasta
                sequence = ""
            
            else:
                sequence = sequence + line.strip()

        if seqname in chromosomes:
            raise Exception(f"Duplicate sequence name {seqname} found in {fasta_file}")

        chromosomes[seqname] = sequence
        debug(f"Added {seqname} {len(sequence)} bp")

    return chromosomes
        
        

def read_gtf(gtf_file, gene_filter):
    debug(f"Reading GTF {gtf_file} with genes {gene_filter}")

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

            # For the gene name and gene id we need to delve into the
            # extended comments
            comments = sections[8].split(";")
            gene_id=None
            gene_name=None
            transcript_id=None
            transcript_name=None

            for comment in comments:
                if comment.strip().startswith("gene_id"):
                    gene_id=comment[8:].replace('"','').strip()
                
                if comment.strip().startswith("gene_name"):
                    gene_name=comment.strip()[10:].replace('"','').strip()
                    
                if comment.strip().startswith("transcript_id"):
                    transcript_id=comment[15:].replace('"','').strip()
                          
                if comment.strip().startswith("transcript_name"):
                    transcript_name=comment.strip()[17:].replace('"','').strip()
                    
            if gene_id is None and gene_name is None:
                warn(f"No gene name or id found for exon at {chrom}:{start}-{end}")
                continue

            if transcript_id is None and transcript_name is None:
                warn(f"No transcript name or id found for exon at {chrom}:{start}-{end}")
                continue


            if gene_id is None:
                gene_id = gene_name

            if gene_name is None:
                gene_name = gene_id

            if gene_filter is not None:
                if not (gene_name == gene_filter or gene_id == gene_filter) :
                    continue
            
            if transcript_id is None and transcript_name is None:
                warn(f"No name or id found for transcript at {chrom}:{start}-{end}")
                continue

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




# To convert the start/stop exon locations to a splice site tuple which looks like
# ((ex1_end),(ex2_start,ex2_end),(ex3_start))
def convert_splice_pattern(exon_list):
    sorted_exons = sorted(exon_list)
    splice_pattern = []
    
    for i,exon in enumerate(sorted_exons):
        if i==0:
            splice_pattern.append((exon[1],)) # Just the end of the first exon

        elif i==len(sorted_exons)-1:
            splice_pattern.append((exon[0],)) # Just the start of the last exon

        else:
            splice_pattern.append((exon[0],exon[1])) # Both ends of the internal exons

    return tuple(splice_pattern)


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
        "--outfile","-o",
        help="The file to write the output count table to",
        default="nexons_output.txt"
    )

    parser.add_argument(
        "--both_out","-b",
        help="Write out to gtf and custom nexons format",
        action="store_true"
    )

    parser.add_argument(
        "--gtf_out","-t",
        help="Write to gtf file instead of custom nexons format",
        action="store_true"
    )

    parser.add_argument(
        "--flexibility","-f",
        help="How many bases different can exon boundaries be and still merge them",
        default=5, 
        type=int
    )

    parser.add_argument(
        "--mincount","-m",
        help="What is the minimum number of observations for any given variant to report",
        default=2, 
        type=int
    )

    parser.add_argument(
        "--mincoverage","-c",
        help="What is the minimum proportion of the gene which must be covered by a variant",
        default=0.1,
        type=float
    )

    parser.add_argument(
        "--mapthreshold","-mt",
        help="What is the minimum proportion of the sequence that must be mapped to a transcript",
        default=0.1,
        type=float
    )

    parser.add_argument(
        "--minexons","-e",
        help="What is the minimum number of exons which must be present to count a transcript",
        default=2,
        type=int
    )

    parser.add_argument(
        "--gene","-g",
        help="The name or ID of a single gene to quantitate"
    )

    parser.add_argument(
        "--direction","-d",
        help="The directionality of the library (none, same, opposing)",
        default="none"
    )

    parser.add_argument(
        "--no_chr_prefix","-nc",
        help="Don't add a 'chr' prefix in temp bed file (required for SIRV sequences where chr doesn't start with 'chr'",
        action="store_true"
    )

    parser.add_argument(
        "--verbose","-v",
        help="Produce more verbose output",
        action="store_true"
    )

    parser.add_argument(
        "--verbose_proportions","-vp",
        help="Print out all proportions of sequences mapped",
        action="store_true"
    )

    parser.add_argument(
        "--report_all",
        help="Report all transcripts, including seeded ones with no counts. Only works with --gtf_out or --both_out. Writes to a separate file named match_info.txt",
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

