
'''
    This script produces a csv file of mapped read counts from a sam file.

    To use:
    python counts.py infile.sam

    out put:
    infile_counts.csv

    out file format:
    gene1,# of mapped reads
    gene2,# of mapped reads
'''

import sys

sam_files = sys.argv[1:]
gene_dic = {}

print "reading in sam header"
with open(sam_files[0], 'r') as in_1:
    for line in in_1:
        if "HWI-" in line:
            break
        else:
            gene_dic[line.split('\t')[1].strip('SN:')]=0
print "done"

bad_seq = "N" * 100

for file in sam_files:
    print "processing file "+file
    count = 0
    with open(file,'r') as infile:
        for line in infile:
            if 'SN:' in line:
                continue
            else:
                elements = line.split('\t')
                count +=1
                if count%1000 == 0:
                    print "on line "+str(count)
                if bad_seq in elements:
                    continue
                else:
                    try:
                        gene_dic[elements[2]]+=1 
                    except KeyError:
                        pass
    with open(file.strip('.sam')+'_counts.csv','w') as outfile:
        print "writing out file"
        outfile.write('%s,%s\n' % ('contig',file))
        for j, contig in enumerate(gene_dic):
            outfile.write('%s,%s\n' % (contig,gene_dic[contig]))

