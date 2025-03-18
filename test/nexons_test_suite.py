#!/usr/bin/env python3

# This script has a series of test cases to validate the basic functionality
# of the main nexons.py script.
# 
# It has tests for the various sub-functions to ensure they are handling
# both main and corner cases correctly.

# We need to add the folder containing nexons.py to the path so we can
# import it

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent.parent))

from nexons import get_possible_genes, match_exons, gene_matches,read_gtf, build_index


def main():
    print("TESTING EXON MATCHING\n---------------------")
    test_exon_matching()

    print("\nTESTING GENE MATCHING\n---------------------")
    test_gene_matching()

    print("\nTESTING FEATURE RETRIEVAL\n-------------------------")
    test_feature_retrieval()

    print("\nTESTING PERCENTILE MATCHING\n-----------------------------")
    test_percentile_matching()



def failed(message):
    print('\033[91m'+"FAILED: "+'\033[0m'+message, file=sys.stderr)

def passed(message):
    print('\033[92m'+"PASSED: "+'\033[0m'+message)


def test_gene_matching():
    gene = {'name': 'TEST', 'id': 'ENSG00000000001', 'chrom': '1', 'start': 10, 'end': 800, 'strand': '+', 'transcripts': {'ENST00000000001': {'name': 'TEST-101', 'id': 'ENST00000000001', 'chrom': '1', 'start': 10, 'end': 800, 'strand': '+', 'exons': [[10,100], [600,800]]}, 'ENST00000000002': {'name': 'TEST-102', 'id': 'ENST00000000002', 'chrom': '1', 'start': 10, 'end': 500, 'strand': '+', 'exons': [[10, 100], [200,500]]}}}

    # Unique hit
    answer = gene_matches([[10,100],[200,500]],gene,0,0)
    if not answer[0] is not None:
        failed("Unique hit not reporting transcript")
    elif not answer[0] == "ENST00000000002":
        failed("Unique hit not reporting correct transcript")
    elif answer[1] != "unique":
        failed("Unique hit not repoted as unique")
    else:
        passed("Unique gene match OK")


    # Partial hit
    answer = gene_matches([[10,100],[200,300]],gene,0,0)
    if not answer[0] is not None:
        failed("Partial hit not reporting transcript")
    elif not answer[0] == "ENST00000000002":
        failed("Partial hit not reporting correct transcript")
    elif answer[1] != "partial":
        failed("Partial hit not repoted as partial")
    else:
        passed("Partial gene match OK")


    # Multi hit - should match multiple transcripts
    answer = gene_matches([[10,100]],gene,0,0)
    if not answer[0] is not None:
        failed("Multi hit not reporting transcript")
    elif answer[1] != "multi":
        failed("Multi hit not repoted as multi")
    else:
        passed("Multi transcript match OK")


    # Intron match
    answer = gene_matches([[500,600]],gene,0,0)
    if answer[0] is not None:
        failed("Intron matching incorrectly reporting transcript")
    elif answer[1] != "intron":
        failed("Intron match not repoted as intron")
    else:
        passed("Intron match OK")




def test_exon_matching():

    # match_exons( exons, transcript, flex, endflex)

    # Perfect match
    answer = match_exons([[10,20],[30,40],[50,60]],[[10,20],[30,40],[50,60]],0,0)
    if not answer[0]:
        failed("Perfect match reported as not matching")

    elif answer[1]:
        failed("Perfect match reported as partial")

    elif any([x!=0 for x in answer[2]]):
        failed("Nonzero endflex for perfect match")

    elif any([x!=0 for x in answer[3]]):
        failed("Nonzero innerflex for perfect match")

    else:
        passed("Perfect match OK")

    # Match with internal flex
    answer = match_exons([[10,25],[30,40],[45,60]],[[10,20],[30,40],[50,60]],5,0)
    if not answer[0]:
        failed("Internal flex reported as not matching")

    elif answer[1]:
        failed("Perfect internal flex reported as partial")

    elif any([x!=0 for x in answer[2]]):
        failed("Nonzero endflex for internal flex")

    elif answer[3] != [5,0,0,-5]:
        failed("Incorreect flex values for internal flex")

    else:
        passed("Internal flex OK")


    # Match with too much internal flex
    answer = match_exons([[10,25],[30,40],[45,60]],[[10,20],[30,40],[50,60]],4,0)
    if answer[0]:
        failed("Too much internal flex reported as matching")

    else:
        passed("Too much internal flex OK")


    # Match with end flex
    answer = match_exons([[5,20],[30,40],[50,65]],[[10,20],[30,40],[50,60]],0,5)
    if not answer[0]:
        failed("End flex reported as not matching")

    elif answer[1]:
        failed("Perfect end flex reported as partial")

    elif answer[2] != [-5,5]:
        failed("Incorreect endflex values for end flex")

    elif any([x!=0 for x in answer[3]]):
        failed("Nonzero innerflex for end flex")

    else:
        passed("End flex OK")


    # Match with too much end flex
    answer = match_exons([[5,20],[30,40],[50,65]],[[10,20],[30,40],[50,60]],0,4)
    if answer[0]:
        failed("Too much end flex reported as matching")

    else:
        passed("Too much end flex OK")


    # Partial match
    answer = match_exons([[35,40],[50,60]],[[10,20],[30,40],[50,60]],0,0)
    if not answer[0]:
        failed("Partial match reported as not matching")

    elif not answer[1]:
        failed("Partial match reported as perfect")

    elif any([x!=0 for x in answer[2]]):
        failed("Nonzero endflex for partial match")

    elif any([x!=0 for x in answer[3]]):
        failed("Nonzero innerflex for partial match")

    else:
        passed("Partial match OK")


    # Internal exon match
    answer = match_exons([[32,38]],[[10,20],[30,40],[50,60]],0,0)
    if not answer[0]:
        failed("Internal exon match reported as not matching")

    elif not answer[1]:
        failed("Internal exon match reported as perfect")

    elif answer[2]:
        failed("Endflex reported for internal exon match")

    elif answer[3]:
        failed("Innerflex reported for internal exon match")

    else:
        passed("Internal exon match OK")


    # Intron match
    answer = match_exons([[25,28]],[[10,20],[30,40],[50,60]],0,0)
    if answer[0]:
        failed("Intron match reported as properly matching")

    elif answer[2]:
        failed("Endflex reported for intron match")

    elif answer[3]:
        failed("Innerflex reported for intron match")

    else:
        passed("Intron match OK")


    # Exon to intron match
    answer = match_exons([[35,45]],[[10,20],[30,40],[50,60]],0,0)
    if answer[0]:
        failed("Exon to Intron match reported as properly matching")

    elif answer[2]:
        failed("Endflex reported for exon to intron match")

    elif answer[3]:
        failed("Innerflex reported for exon to intron match")

    else:
        passed("Exon to Intron match OK")


    # Mismatch
    answer = match_exons([[100,150]],[[10,20],[30,40],[50,60]],0,0)
    if answer[0]:
        failed("Mismatch reported as matching")
        
    else:
        passed("Mismatch exons OK")
        

def test_feature_retrieval():
    # test.gtf is human GRCh38v113 chr3 between 63723373-64479723

    # TSL Filtering
    all_genes_tsl1 = read_gtf(Path(__file__).parent/"test.gtf",1)
    if not len(all_genes_tsl1["ENSG00000163635"]["transcripts"].keys()) == 8:
        failed("ENSG00000163635 at TSL1 didn't have 8 transcripts")
    else:
        passed("TSL1 Filtering OK")

    all_genes_tsl5 = read_gtf(Path(__file__).parent/"test.gtf",5)    
    if not len(all_genes_tsl5["ENSG00000163635"]["transcripts"].keys()) == 13:
        failed("ENSG00000163635 at TSL5 didn't have 13 transcripts")
    else:
        passed("TSL5 Filtering OK")

    index = build_index(all_genes_tsl1)
    index5 = build_index(all_genes_tsl5)

    # Single fetch
    answer = get_possible_genes(index,"3",64100000,64200000,"A")
    if not len(answer) == 1:
        failed("No gene found for single fetch")
    elif [x for x in answer][0] != 'ENSG00000163637':
        failed("Wrong gene retrieved for single fetch")
    else:
        passed("Simple gene fetch OK")

    # Directional fetch
    answer = get_possible_genes(index,"3",64100000,64200000,"R")
    if not len(answer) == 1:
        failed("No gene found for directional fetch")
    elif [x for x in answer][0] != 'ENSG00000163637':
        failed("Wrong gene retrieved for directional fetch")
    else:
        passed("Directional gene fetch OK")

    # Wrong Directional fetch
    answer = get_possible_genes(index,"3",64100000,64200000,"F")
    if not len(answer) == 0:
        failed("Incorrect hit found for wrong directional fetch")
    else:
        passed("Wrong directional gene fetch OK")

    # Multiple Fetch
    answer = get_possible_genes(index5,"3",64190000,64190001,"A")
    if not len(answer) == 2:
        failed("Didn't find two hits for multiple fetch")
    else:
        passed("Multiple fetch OK")
    
    # Exact start
    answer = get_possible_genes(index,"3",63819299,63819299,"F")
    if not len(answer) == 1:
        failed("No gene found for exact start fetch")
    else:
        passed("Exact start fetch OK")

    # Exact end
    answer = get_possible_genes(index,"3",64200965,64200965,"F")
    if not len(answer) == 1:
        failed("No gene found for exact end fetch")
    else:
        passed("Exact end fetch OK")


def test_percentile_matching():
    # Full start to end match
    answer = match_exons([[100,200],[300,400],[500,600]],[[100,200],[300,400],[500,600]],0,0)
    if not answer[4] == 0:
        failed("Incorrect start for full match")
    elif not answer[5] == 100:
        failed("Incorrect end for full match")
    else:
        passed("Percentile full match OK")
    

    # Start to mid match
    answer = match_exons([[100,200],[300,350]],[[100,200],[300,400],[500,600]],0,0)
    if not answer[4] == 0:
        failed("Incorrect start for start mid match")
    elif abs(50-answer[5])>1:
        failed("Incorrect end for start mid match")
    else:
        passed("Percentile start mid match OK")


    # Mid to end match
    answer = match_exons([[350,400],[500,600]],[[100,200],[300,400],[500,600]],0,0)
    if not answer[5] == 100:
        failed("Incorrect end for mid end match")
    elif abs(50-answer[4])>1:
        failed("Incorrect start for start mid match")
    else:
        passed("Percentile mid end match OK")


    # Inner match
    answer = match_exons([[350,400],[500,550]],[[100,200],[300,400],[500,600]],0,0)
    if abs(50-answer[4])>1:
        failed("Incorrect start for inner match")

    if abs(83-answer[5])>1:
        failed("Incorrect end for inner match")

    else:
        passed("Percentile Inner match OK")


    # Single exon match
    answer = match_exons([[550,600]],[[100,200],[300,400],[500,600]],0,0)
    if abs(answer[4]-83) > 1:
        failed("Incorrect start single exon match")
    elif not answer[5] == 100:
        failed("Incorrect end for single exon match")
    else:
        passed("Percentile single exon match OK")


    # Single exon match 2
    answer = match_exons([[100,150]],[[100,200],[300,400],[500,600]],0,0)
    if abs(answer[5]-17) > 1:
        failed("Incorrect end single exon match2")
    elif not answer[4] == 0:
        failed("Incorrect start for single exon match2")
    else:
        passed("Percentile single exon match2 OK")


    # Middle exon match
    answer = match_exons([[300,400]],[[100,200],[300,400],[500,600]],0,0)
    if abs(answer[4]-33) > 1:
        failed("Incorrect middle exon start match")
    elif abs(answer[5]-67) > 1:
        failed("Incorrect middle exon end match")
    else:
        passed("Percentile middle exon match OK")


    # Single exon transcript
    answer = match_exons([[100,150]],[[50,200]],0,0)
    if abs(answer[4]-33) > 1:
        failed("Incorrect start for single exon transcript")
    elif abs(answer[5]-67) > 1:
        failed("Incorrect end for single exon transcript")
    else:
        passed("Percentile single exon transcript OK")
    

if __name__ == "__main__":
    main()