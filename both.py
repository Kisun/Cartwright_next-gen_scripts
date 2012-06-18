#!/usr/bin/python

'''
	both.py processes two fastq files generated using Illumina's pipeline CASAVA 1.8 and higher, 
	one containing R1 reads and another containing R2 reads and returns two files containg only reads 
	that are present in both of the input files.

	To use:

	python both.py R1.fastq R2.fastq 

	out put:
	R1.fastq.both 
	R2.fastq.both

	This script uses the screed module: https://github.com/ctb/screed
'''

import screed
import sys

R1_IN = sys.argv[1]
R2_IN = sys.argv[2]

screed.read_fastq_sequences(R1_IN)
screed.read_fastq_sequences(R2_IN)

DB_R1 = screed.ScreedDB(R1_IN+'_screed')
DB_R2 = screed.ScreedDB(R2_IN+'_screed')

with open(R1_IN+'.both','w') as R1_OUT:
    with open(R2_IN+'.both','w') as R2_OUT:
        for record, thing in DB_R1.iteritems():
            try:
                match = DB_R2[thing['name'].replace(" 1:"," 2:")]
            except KeyError:
                continue
            R1_OUT.write('@%s %s\n%s\n+\n%s\n' % (thing['name'],thing['annotations'],thing['sequence'],thing['accuracy']))      
            R2_OUT.write('@%s %s\n%s\n+\n%s\n' % (match['name'],match['annotations'],match['sequence'],match['accuracy']))
