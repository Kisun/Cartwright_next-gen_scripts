'''


Calculate percent diferences in up/down regulation for mapped read counts. 

Mapped read counts are normalized by dividing the raw count by the gene 
length multiplied by the total number of mapped reads. Percent differences 
are computed as the absolute value of the difference between normalized 
mapped read counts, divided by the larger mapped read count, multiplied 
by 100. 

usage: python percent.py filein percent difference normalize 

filein: a .csv file of mapped read counts --> "gene, counts sample1, counts
sample2, gene length"
percent: This script will out put a file with all normalized counts and 
percentages. A percent cut off should be specified in this argument to produce
an additional file with only genes whose percent differences are equal to and
above the percent cut off. 
difference: a difference cut off between normalized mapped read counts.
normalize: specify 'y' if normalization is required, otherwise the mapped read
counts in the input file will be treated as normalized. 



'''
import sys
import csv
from math import fabs as abs, floor
import os

file = csv.reader(open(sys.argv[1]))
percent = int(sys.argv[2])
diff = int(sys.argv[3])
norm=sys.argv[4]
dist = {'10:20':0,'20:30':0,'30:40':0,'40:50':0,'50:60':0,'60:70':0,'70:80':0,'80:90':0,'90:101':0}

gene, sample1, sample2, length  = zip(*file)

if norm == 'y':
    sample1_len = sum([int(i) for i in sample1[1:]])
    sample2_len = sum([int(i) for i in sample2[1:]])

    sample1_norm = [str(int(i[0])/(float(int(i[1])*sample1_len))) for i in zip(sample1[1:], length[1:])]
    sample2_norm = [str(float(int(i[0]))/(float(int(i[1]))*sample2_len)) for i in zip(sample2[1:], length[1:])]
else:
    sample1_norm = sample1[1:]
    sample2_norm = sample2[1:]
    
myfile = os.path.split(sys.argv[1])[1].strip(".csv")
outfile = '%s_norm_%s.csv' %(myfile,str(percent))
rawfile = '%s_norm_raw.csv' %myfile
sfile = '%s_stats.csv'  %myfile

with open(rawfile,'w') as out:
    with open(outfile, 'w') as out2:
        out.write('gene,length,%(1)s raw,%(2)s raw,%(1)s normalized,%(2)s normalized,difference,percent\n' %{'1':sample1[0], '2':sample2[0]})
        out2.write('gene,length,%(1)s raw,%(2)s raw,%(1)s normalized,%(2)s normalized,difference,percent\n' %{'1':sample1[0], '2':sample2[0]})
        for j in zip(gene[1:],length[1:],sample1[1:],sample2[1:],sample1_norm,sample2_norm):
            try:
                sample1_m_sample2 = (abs(float(j[4]) - float(j[5]))/max(float(j[4]),float(j[5])))*100 
                sample1_sample2_diff = abs(float(j[4]) - float(j[5]))
            except ZeroDivisionError:
                sample1_m_sample2 = 0
                sample1_sample2_diff = abs(float(j[4]) - float(j[5]))
            for i,k in enumerate(dist.keys()):
	        myrange = k.split(':')
	        if int(floor(sample1_m_sample2)) in range(int(myrange[0]),int(myrange[1])):
	            dist[k]+=1
	    if sample1_m_sample2 >= percent:
                if sample1_sample2_diff >= diff:
                    out2.write(('%s,%s,%s' % (','.join(j), str(sample1_sample2_diff), str(sample1_m_sample2))))  
                    out2.write('\n')
            out.write('%s,%s,%s' % (','.join(j), str(sample1_sample2_diff), str(sample1_m_sample2)))
            out.write('\n')
	Totaldist = dict([('>=%s'%m,sum(dist[n] for n in dist.keys() if int(n.split(':')[0])>=m)) for m in range(10,100,10)])
#	print Totaldist

for i, j in enumerate(Totaldist.keys()):
    print '%s\t%s' %(j, Totaldist[j])
