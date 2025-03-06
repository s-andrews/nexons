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

from nexons import get_possible_genes, match_exons


def main():
    test_exon_matching()
    test_feature_retrieval()

def failed(message):
    print("FAIL: "+message, file=sys.stderr)

def passed(message):
    print("PASSED: "+message)


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
        

    pass


def test_feature_retrieval():
    pass


if __name__ == "__main__":
    main()