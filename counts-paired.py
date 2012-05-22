import sys

'''
    This script produces a csv file of mapped read counts from a sam file. 
    Reads are only counted if both pais map to the same scaffold and if both
    reads are sequenced sucessfully (neither pair is composed of all N's).

    To use:
    python counts-paired.py infile.sam

    out put:
    infile_counts.csv

    out file format:
    gene1,# of mapped reads
    gene2,# of mapped reads

'''

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
        AT1 = infile.readline()
        line1 = infile.readline()
        line2 = infile.readline()
        while line2:
            line1 = infile.readline()
            line2 = infile.readline()
            if 'SN:' in line1 and 'SN:' in line2 or '@PG' in line2:
                continue
            else:
                elements1 = line1.split('\t')
                elements2 = line2.split('\t')
                count +=2
                if count%1000 == 0:
                    print "on line "+str(count)
                if bad_seq in elements1 or bad_seq in elements2:
                    continue
                else:
                    try:
                        if elements1[0] == elements2[0] and elements1[2] == elements2[2]:
                            try:
                                gene_dic[elements1[2]]+=1 
                            except KeyError:
                                pass
                    except IndexError:
                        pass
    with open(file.strip('.sam')+'_counts.csv','w') as outfile:
        print "writing out file"
        outfile.write('%s,%s\n' % ('contig',file))
        for j, contig in enumerate(gene_dic):
            outfile.write('%s,%s\n' % (contig,gene_dic[contig]))

