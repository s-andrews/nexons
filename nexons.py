#!/usr/bin/env python
import argparse
import warnings
import tempfile
import subprocess
import os

verbose = False
gtf_out = False

def main():
    options = get_options()

    global verbose
    global gtf_out
    verbose = options.verbose
    gtf_out = options.gtf_out
    
    # for now, we're going to read the gtf multiple times, as first we want the genes, then transcripts, then exons. 
    # If the gtf was ordered, we could do it in one pass, but safer to do it 3 times - might be slow though...
    genes = read_gtf_genes(options.gtf, options.gene)   
    
    genes_transcripts = read_gtf_transcripts(options.gtf, options.gene, genes)
    
    genes_transcripts_exons = read_gtf_exons(options.gtf, options.gene, genes_transcripts)


    # now all the exons have been imported from the gtf, we can convert the splice patterns to a compatible format.
    for gene_id in genes_transcripts_exons:
        for transcript_id in genes_transcripts_exons[gene_id]["transcripts"]:
            exon_set = genes_transcripts_exons[gene_id]["transcripts"][transcript_id]["exons"]
            splices = sort_splice_pattern(exon_set)
            genes_transcripts_exons[gene_id]["transcripts"][transcript_id]["splice_patterns"] = splices
     
    if verbose:
        print(f"Found {len(genes.keys())} genes")
       # print(f"Found {len(transcripts.keys())} transcripts")
    
    if len(genes.keys())==0:
        raise Exception("No genes to process")

    chromosomes = read_fasta(options.fasta)
    if verbose:
        print(f"Found {len(chromosomes.keys())} chromosomes")

    if len(chromosomes.keys())==0:
        raise Exception("No chromosomes found")

    # pass in splice_patterns as an argument to process_bam_file

    quantitations = {}
    for bam_file in options.bam:
        quantitations[bam_file] = process_bam_file(genes, chromosomes, bam_file, options.minexons, options.mincoverage)

    collated_splices = collate_splice_variants(quantitations,options.flexibility, genes_transcripts_exons)
    quantitations = collated_splices[0]
    splice_info = collated_splices[1]

    if gtf_out:
        if options.outfile == "nexons_output.txt":
            outfile = "nexons_output.gtf"
        else:
            outfile = options.outfile
        write_gtf_output(quantitations, genes, outfile, options.mincount, splice_info)
    else:
        write_output(quantitations, genes, options.outfile, options.mincount, splice_info)


def collate_splice_variants(data, flexibility, genes_transcripts_exons):
    # The structure for the data is 
    # data[bam_file_name][gene_id][splicing_structure] = count
    # 
    # We will have all genes in all BAM files, but might
    # not have all splice forms in all BAM files
    # 
    # Here we aim to produce a reduced set of splice variants
    # by combining variants which differ by only a few positions

    if verbose:
        print("Merging similar variants")

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

        unknown_transcript = 1
        #splice_counts = {} 
       
        splice_counts[gene] = {}
       
       # simplifying the splice counts - genes_transcripts_exons[gene]["transcripts"] contains a load of info for that transcript - a dict of start, stop, exons, splice_patterns (maybe it doesn't need that but we'll leave it for now).
        for transcript_id in genes_transcripts_exons[gene]["transcripts"]:

            splice_counts[gene][genes_transcripts_exons[gene]["transcripts"][transcript_id]["splice_patterns"]] = {
                "transcript_id": transcript_id,
                "count": 0
            }
       
# get all the splice_patterns

        for bam in bam_files:
            these_splices = data[bam][gene].keys()

            for splice in these_splices:
                #if not splice in splice_counts:
                if not splice in splice_counts.keys():
                    splice_counts[gene][splice] = {
                        "transcript_id": unknown_transcript,
                        "count": 0
                    }
                    unknown_transcript += 1
                
                splice_counts[gene][splice]["count"] += data[bam][gene][splice]
        # now this is in a slightly different format for accessing the counts
        
        # Now we can collate the splice counts.  We'll get
        # back a hash of the original name and the name we're
        # going to use.
        # 
        # We'll give the function the ordered set of splices so
        # it always uses the most frequently observed one if 
        # there is duplication - 
        # but we want the known transcripts at the start.... come back to this and order by most frequent TODO!!
        
        splice_name_map = create_splice_name_map(splice_counts[gene].keys(),flexibility) # this now has 
        #splice_name_map = create_splice_name_map(sorted(splice_counts.keys(), key=lambda x:splice_counts[x], reverse=True),flexibility)


        # From this map we can now build up a new merged set
        # of quantitations which put together the similar splices
        for bam in bam_files:
            these_splices = data[bam][gene].keys()

            for splice in these_splices:
                used_splice = splice_name_map[splice]
                if not used_splice in merged_data[bam][gene]:
                    merged_data[bam][gene][used_splice] = 0
                
                merged_data[bam][gene][used_splice] += data[bam][gene][splice]

    return [merged_data, splice_counts]



def create_splice_name_map(splices, flexibility):
    if verbose:
        print(f"Merging {len(splices)} different splice sites")
    # This takes an ordered list of splice strings and matches
    # them on the basis of how similar they are.  We work 
    # our way down the list trying to match to previous strings

    # This is what we'll give back.  Every string will be in this
    # and we'll match it either to itself or a more popular 
    # equivalent string

    map_to_return = {}

    # We'll store the good split strings in here so they're easier
    # to compare to
    sized_segments = {}

    for splice in splices:
        # Split it into segments and then parse the numbers out of these
        segments = splice.split(":")

        parsed_segments = []

        for segment in segments:
            parsed_segments.append([int(x) for x in segment.split("-")])

        # Now we try to map the parsed segment to existing segments

        # If there are no splices with this size then we accept it and
        # move on
        if len(parsed_segments) not in sized_segments:
            sized_segments[len(parsed_segments)] = [{"segments": parsed_segments, "string":splice}]
            map_to_return[splice] = splice
            continue

        # We now test the existing parsed segments to see if they're
        # close enough to this one.
        found_hit = False
        for test_segment in sized_segments[len(parsed_segments)]:
            too_far = False
            # Iterate though the segments
            for i in range(len(parsed_segments)):
                # Iterate through the start/end positions
                for j in range(len(parsed_segments[i])):
                    if abs(parsed_segments[i][j]-test_segment["segments"][i][j]) > flexibility:
                        too_far = True
                        break
                if too_far:
                    break
            if not too_far:
                # We can use this as a match
                #print(f"Merged:\n{splice}\ninto\n{test_segment['string']}\n\n")
                map_to_return[splice] = test_segment["string"]
                found_hit = True
                break
        
        if not found_hit:
            # Nothing was close enough, so enter this as a new reference
            sized_segments[len(parsed_segments)].append({"segments": parsed_segments, "string":splice})
            map_to_return[splice] = splice

    if verbose:
        dedup_splices = set()
        for x in map_to_return.values():
            dedup_splices.add(x)
        
        print(f"Produced {len(dedup_splices)} deduplicated splices")
    
    return map_to_return



def write_output(data, gene_annotations, file, mincount, splice_info):
    # The structure for the data is 
    # data[bam_file_name][gene_id][splicing_structure] = count
    # 
    # We will have all genes in all BAM files, but might
    # not have all splice forms in all BAM files

    bam_files = list(data.keys())
 
    with open(file,"w") as outfile:
        # Write the header
        header = ["Gene ID", "Gene Name","Chr","Strand","SplicePattern", "Transcript id"]
        header.extend(bam_files)
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

            # Now we can go through the splices for all BAM files
            for splice in splices:
                line_values = [gene,gene_annotations[gene]["name"],gene_annotations[gene]["chrom"],gene_annotations[gene]["strand"],splice]

                line_above_min = False
                for bam in bam_files:
                    if splice in data[bam][gene]:
                        
                        line_values.append(splice_info[gene][splice]["transcript_id"])
                        line_values.append(data[bam][gene][splice]) # adding the count
                       # line_values.append(splice_info[gene][splice]["count"]) # this count doesn't work
                        if data[bam][gene][splice] >= mincount:
                            line_above_min = True
                    
                    else:
                        line_values.append(0)

                if line_above_min:
                    lines_written += 1
                    outfile.write("\t".join([str(x) for x in line_values]))
                    outfile.write("\n")

    if verbose:
        print(f"Wrote {lines_written} splices to {file}")


def write_gtf_output(data, gene_annotations, file, mincount, splice_info):
    # The structure for the data is 
    # data[bam_file_name][gene_id][splicing_structure] = count
    # 
    # We will have all genes in all BAM files, but might
    # not have all splice forms in all BAM files

    bam_files = list(data.keys())
 
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

                line_values = [gene_annotations[gene]["chrom"], "nexons", "transcript", splice_start, splice_end, ".", gene_annotations[gene]["strand"], 0]

                splice_text = "splicePattern " + splice

                line_above_min = False
                for bam in bam_files:
                    if splice in data[bam][gene]:
                        
                        transcript_text = "transcript_id " + str(splice_info[gene][splice]["transcript_id"])
                        
                        attribute_field = gtf_gene_text + "; " + transcript_text + "; " + splice_text
                        
                        line_values[5] = data[bam][gene][splice] # adding the count                        
                        line_values.append(attribute_field)

                        if data[bam][gene][splice] >= mincount:
                            line_above_min = True
                    
                    else:
                        line_values.append(0)

                if line_above_min:
                    lines_written += 1
                    outfile.write("\t".join([str(x) for x in line_values]))
                    outfile.write("\n")

    if verbose:
        print(f"Wrote {lines_written} splices to {file}")


def process_bam_file(genes, chromosomes, bam_file, min_exons, min_coverage):
    counts = {}

    # we need a set of splices that are passed in - that need to be in a dictionary of gene_ids as there may be 
    # multiple genes we're looking for.


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

        fasta_file = tempfile.mkstemp(suffix=".fa", dir="/dev/shm")
        sequence = chromosomes[gene["chrom"]][(gene["start"]-1):gene["end"]]

        with open(fasta_file[1],"w") as out:
            out.write(f">{gene['name']}\n{sequence}\n")

 
        gene_counts = {}# this is where we seed the transcripts I think 
        reads = get_reads(gene,bam_file)

        if verbose:
            print(f"Got {len(reads)} reads for {gene['name']} from {bam_file}")

        for read_id in reads.keys():

            # See if we need to print out a progress message
            progress_counter += 1
            if verbose and progress_counter % 100 == 0:
                print("Processed "+str(progress_counter)+" reads  from "+bam_file)
                ## FOR TESTING ONLY ###
                #break

            
            # Running chexons will fail if there are no matches between 
            # the read and the gene (or too short to make a match segment)
            # We therefore need to catch this and skip any reads which 
            # trigger it

            try:
                segment_string = get_chexons_segment_string(reads[read_id],fasta_file[1],gene, min_exons, min_coverage)

                # This might be None if the transcript didn't pass the filters
                # we've put in place
                if segment_string is not None:
                    if not segment_string in gene_counts:
                        gene_counts[segment_string] = 0
                
                    gene_counts[segment_string] += 1
                            
            except Exception as e:
                print(f"[WARNING] Chexons failed with {e}")


        # Clean up the gene sequence file
        os.unlink(fasta_file[1])
        
        counts[gene["id"]] = gene_counts

    return counts

def get_chexons_segment_string (sequence, genomic_file, gene, min_exons, min_coverage):

    position_offset = gene["start"]
    # We need to write the read into a file
    read_file = tempfile.mkstemp(suffix=".fa", dir="/dev/shm")

    with os.fdopen(read_file[0],"w") as out:
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
            start = int(start_end[0])+position_offset-1
            end = int(start_end[-1])+position_offset-1

            locations.append([start,end])

    # Check that we've got enough exons to keep this
    if len(locations) < min_exons:
        return None

    # Check that we've got enough coverage
    if abs(locations[0][0] - locations[-1][-1]) / abs(gene["end"]-gene["start"]) < min_coverage:
        return None 

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
    os.remove(read_file[1])


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
    
    bed_file = tempfile.mkstemp(suffix=".bed", dir="/dev/shm")

    with open(bed_file[1],"w") as out:
        ## TODO: work out how to handle chr names (chr prefix)
        #out.write(f"chr{gene['chrom']}\t{gene['start']}\t{gene['end']}\n")
        out.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\n") #Laura 

    # Now we can run samtools to get the data.  We require that the
    # read overlaps the relevant region, but we also require that
    # the reads is on the opposing strand to the gene (since that's how
    # the directionality of nanopore reads works)
    #
    # The filter for forward strand reads is -f 16 and reverse is -F 16
    # the library is opposing strand specific so if we have a forward
    # strand gene we want -f and we'd use -F for reverse strand.

    strand_filter_string = "-f" if gene["strand"] == "+" else "-F"

    samtools_process = subprocess.Popen(["samtools","view",bam_file,"-L",bed_file[1],strand_filter_string,"16"], stdout=subprocess.PIPE)
    
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
                seqname = seqname.rstrip() # Laura - removing new line character in fasta
                sequence = ""
            
            else:
                sequence = sequence + line.strip()

        if seqname in chromosomes:
            raise Exception(f"Duplicate sequence name {seqname} found in {fasta_file}")

        chromosomes[seqname] = sequence
        if verbose:
            print(f"Added {seqname} {len(sequence)} bp")

    return chromosomes
        
        



def read_gtf_exons(gtf_file, gene_filter, genes_transcripts):

    if verbose:
        print(f"Reading exons from {gtf_file}")

    with open(gtf_file) as file:
        for line in file:
            # Skip comments
            if line.startswith("#"):
                continue

            sections = line.split("\t")

            if len(sections) < 7:
                warnings.warn("Not enough data from line "+line)
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
                warnings.warn(f"No gene name or id found for exon at {chrom}:{start}-{end}")
                continue

            if gene_id is None:
                gene_id = gene_name

            if gene_name is None:
                gene_name = gene_id

            if gene_filter is not None:
                if not (gene_name == gene_filter or gene_id == gene_filter) :
                    continue
            
            if transcript_id is None and transcript_name is None:
                warnings.warn(f"No name or id found for transcript at {chrom}:{start}-{end}")
                continue

            if transcript_id is None:
                transcript_id = transcript_name

            if transcript_name is None:
                transcript_name = transcript_id

            exons = [start, end]

            genes_transcripts[gene_id]["transcripts"][transcript_id]["exons"].append(exons)
                
    return genes_transcripts

def read_gtf_transcripts(gtf_file, gene_filter, genes):

    if verbose:
        print(f"Reading genes and transcripts from {gtf_file}")

    genes_transcripts = genes

    with open(gtf_file) as file:
        for line in file:
            # Skip comments
            if line.startswith("#"):
                continue

            sections = line.split("\t")

            if len(sections) < 7:
                warnings.warn("Not enough data from line "+line)
                continue

            #if sections[2] != "gene" and sections[2] != "transcript":
             #   continue

            if sections[2] != "transcript":
                continue

            # we can pull out the main information easily enough
            chrom = sections[0]
            start = int(sections[3])
            end = int(sections[4])
            strand = sections[6]
            type = sections[2]

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

            # here we need to check whether gene_name or id is what we're after
            if gene_filter is not None:
                if not (gene_name == gene_filter or gene_id == gene_filter) :
                    continue
                
            transcript_id=None
            transcript_name=None
            for comment in comments:
                if comment.strip().startswith("transcript_id"):
                    transcript_id=comment[15:].replace('"','').strip()
                    
                if comment.strip().startswith("transcript_name"):
                    transcript_name=comment.strip()[17:].replace('"','').strip()
                    
            if transcript_id is None and transcript_name is None:
                warnings.warn(f"No name or id found for transcript at {chrom}:{start}-{end}")
                continue

            if transcript_id is None:
                warnings.warn(f"Using transcript name {transcript_name} as ID")
                transcript_id = transcript_name

            if transcript_name is None:
                warnings.warn(f"Using transcript ID {transcript_id} as name")
                transcript_name = transcript_id
                
            genes_transcripts[gene_id]["transcripts"][transcript_id] = {
                "name": transcript_name,
                "id": transcript_id,
                "chrom": chrom,
                "start": start,
                "end" : end,
                "strand" : strand,
                "exons": []
            }
            
    return genes_transcripts

def read_gtf_genes(gtf_file, gene_filter):

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
                "strand" : strand,
                "transcripts": {}
            }

    return genes

# To convert the start/stop exon locations to compatible splice pattern,
# e.g. 20:35-40:55-65:100 
#def sort_splice_pattern(transcript_dict):
def sort_splice_pattern(exon_list):
    sorted_exons = sorted(exon_list)  
    #sorted_exons = sorted(transcript_dict["exons"])   
    flat_list = [str(sorted_exons[0][1])] # we only want the end of the first exon
    sorted_exons.pop(0)

    if sorted_exons: # check it's not empty

        last_element = str(sorted_exons.pop()[0]) # only want the start of the last exon

        for exon in sorted_exons:
            flat_list.append("-".join(map(str, exon)))

        flat_list.append(last_element)
        
    splice_pattern = ':'.join(flat_list)   
    return splice_pattern


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
        "--verbose","-v",
        help="Produce more verbose output",
        action="store_true"
    )

    return(parser.parse_args())



if __name__ == "__main__":
    main()

