#!/usr/bin/env python
import getopt
import sys
import os
import argparse

input_file = ""

parser = argparse.ArgumentParser(description = '''Nexons has an option to create gtf output files, but these do not separate out the exons on to individual rows. 
This script should take a gtf output file from nexons and create a new gtf file with each exon on a separate line.
This should be compatible with IGV. Information from the first and last splice sites currently get thrown away.''')
parser.add_argument('input_file', type=str, default="", help='nexons output gtf file')
args=parser.parse_args()

input_file = args.input_file

outfile = input_file.replace(".gtf", "_exons_extracted.gtf")
print('\n-----------------------------------------')
print(f'writing to file {outfile}')

line_count = 0

out = open(outfile, "w")

with open(input_file) as f:

    line = f.readline()
    
    while line:   
        line_count = line_count+1    
        out.write(f'{line}')
        
        row = line.split("\t")
    
        info = row[8]
        info_split = info.split(";")
        splices = info_split[2]
        splices_split = splices.split(":")
        splices_split.pop(0) 
        
        for exon in splices_split:
            exon_start_end = exon.split("-")
            if len(exon_start_end) == 2:
            
                out.write(f'{row[0]}\t{row[1]}\texon\t{exon_start_end[0]}\t{exon_start_end[1]}\t{row[5]}\t{row[6]}\t{row[7]}\t{row[8]}')

        try:
            line = f.readline()
           
        except OSError:
            print(f'something went wrong on line {line_count}')
        

out.close()
