'''

This script extracts annotations for OrthoMCL classified orthologs using
OrthoMCL annotation databases (defined in the Databases list).

The input is an OrthoMCL output file.

'''

import sys
import os

mclFile = sys.argv[1]

#indicate databases here
Databases = ['/volumes/nextgen_ftp/mariya/orthomcl/aa_deflines_orthomcl-5.txt','/Volumes/Nextgen_FTP/mariya/orthomcl/aa_deflines_OrthoMCL-3.txt', '/Volumes/Nextgen_FTP/mariya/orthomcl/aa_deflines_OrthoMCL-4.txt', '/Volumes/Nextgen_FTP/mariya/orthomcl/BAE_geneid_anno' ]

with open(mclFile,'r') as infile:
    for line in infile:
        myline = line.split()
        for base in Databases:
            p = os.popen('grep "%s" %s 2>err' % (myline[2].strip('\n')+' | '+myline[1],base))
            query = p.readline()
            if query != '':
                querynew = query.strip('>').strip('\n').split('|')[3]
                myline.insert(3, querynew)
                if querynew == ' ':
                    for data in Databases:
                        p = os.popen('grep "%s" %s 2>err' % (myline[1],data))
                        q = p.readlines()
                        if q != []:
                            print '\t'.join(myline)
                            for element in q:
                                mylist = element.strip('\n').split('|')
                                if mylist[3] != ' ' and mylist[3].strip() != "hypothetical protein":
                                    print "\t%s\t%s|%s\t%s" % (mylist[2], mylist[0].strip('>'), mylist[1], mylist[3])
                            break
                        else:
                            continue
                else:
                    print '\t'.join(myline)
                break
            else:
                continue
