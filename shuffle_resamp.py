'''

This script creates subsampled datasets of paired-end Illumina Reads. 

The Screed library is required to execute this script. Screed can be downloaded
from: git@github.com:ctb/screed.git .

usage: python shuffle_resamp.py filein percet number
filein: a fastq file containing both read pairs.
percent: % of dataset to resample.
number: number of resampled datasets to create. 

'''

#!/usr/bin/python
from __future__ import print_function

import screed
import sys
import os
import logging, logging.handlers
from random import randint 
from math import ceil

def resamp(percent, tempfile, mypath, fileout, db, size, numresamp):
    percent, size = float(percent), float(size)
    print('The size of the original population is ',size)   
    NewSize = (percent/100)*size
    myrange = ceil(size/NewSize)
    NumFiles = dict([(x,(x*NewSize-(NewSize-1),x*NewSize)) for x in range(1,int(myrange)+1)])
    print(percent,'% will be resampled from the population. ',2*NewSize,' sequences will be written into ',numresamp,' files.')
    print(NumFiles)
    i = 0
    with open(tempfile,'r') as temp:
            for line in temp:
                i += 1
                myNum = randint(1,size)
                if i%1000 == 0:
                    print('Processing sequence ',i)
                for x in range(1,int(numresamp)+1):
                    myfile = fileout + '_' + str(x)
                    if NumFiles[x][0] <= myNum and NumFiles[x][1] >= myNum:
                        R1 = myfile + 'R1'
                        R2 = myfile + 'R2'
                        write_R1 = open(os.path.join(mypath,R1),'a')
                        write_R2 = open(os.path.join(mypath,R2),'a')
                        myname = line.strip('\n')
                        match2 = db[myname.replace(" 1:"," 2:")]
                        write_R1.write("@"+db[myname]['name'].replace(" ","_")+'/1\n'+db[myname]['sequence']+'\n'+'+'+'\n'+db[myname]['accuracy']+'\n')
                        write_R2.write('@'+match2['name'].replace(" ","_")+'/2\n'+match2['sequence']+'\n'+'+'+'\n'+match2['accuracy']+'\n')       
                        write_R1.close()
                        write_R2.close()
                    else:
                        continue
        
def log(Mypath, ReadNum, Filein, Fileout, Singleton1, Singleton2,n,j,k,l):
    # Make a global logging object.
    x = logging.getLogger("logfun")
    x.setLevel(logging.DEBUG)

    # This handler writes a log file.
    h1 = logging.FileHandler(os.path.join(Mypath,"shuffle_resamp.log"))
    f = logging.Formatter("%(asctime)s %(lineno)d %(message)s")
    h1.setFormatter(f)
    h1.setLevel(logging.DEBUG)
    x.addHandler(h1)
    
    logfun = logging.getLogger("logfun")
    
    logfun.debug("There are "+str(ReadNum)+" sequences in "+Filein)
    logfun.debug(str(n)+" sequences processed.")
    logfun.debug(str(j)+ " sequences written to "+Fileout)
    logfun.debug(str(k)+ " sequences written to R1 "+Singleton1)
    logfun.debug(str(l)+ " sequences written to R2 "+Singleton2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mypath, File = os.path.split(sys.argv[1])
Filein = sys.argv[1]
Fileout = Filein + '.shuffled'
Singleton1 = Filein + '.singletons_R1'
Singleton2 = Filein + '.singletons_R2'
NumResamp = sys.argv[3] #number of resamp datasets to produce

try:
    Percent = sys.argv[2]
    try:
        tempfile = open(os.path.join(Mypath,'temp'),'r')
        tempfile.close()
    except IOError:
        pass
    if not Percent.isdigit():
        print ('You did not enter a valid percent cutoff for resamplaing.')
        sys.exit()
except IndexError:
    Percent = 100
   
screeddb = Filein + '_screed'

if os.path.exists(screeddb):
    db = screed.ScreedDB(screeddb)
    if Percent != 100 and os.path.exists(os.path.join(Mypath,'temp')):
        print("Screed database and shuffled files already exists proceeding to resampling.")
        resamp(Percent, os.path.join(Mypath,'temp'), Mypath, Fileout, db, sys.argv[3], NumResamp)	  
        os.remove(os.path.join(Mypath,'temp'))
        print('Shuffleing done with resampleing. No log file written.')
        sys.exit()
    else:      
        print("Screed database already exists proceeding to shuffleing.")
else:
    print("Creating screed database.")
    db = screed.read_fastq_sequences(Filein)

print(dict(db))
ReadNum = len(db)


N='name'
S='sequence'
A='accuracy'

n, j, k, l = 0, 0, 0, 0

with open(Fileout, 'w') as outfile:
    with open(Singleton1, 'w') as R1:
        with open(Singleton2, 'w') as R2:
            for record, thing in db.iteritems():
                n+=1
		if " 1:" in thing[N]:
                    try:
                        match2 = db[thing[N].replace(" 1:"," 2:")]
                        outfile.write("@"+thing[N].replace(" ","_")+'/1\n'+thing[S]+'\n'+'+'+'\n'+thing[A]+'\n'+'@'+match2[N].replace(" ","_")+'/2\n'+match2[S]+'\n'+'+'+'\n'+match2[A]+'\n') 
                    	j+=2
			if Percent != 100:
                            with open(os.path.join(Mypath,'temp'),'a') as tempfile:
                                tempfile.write(thing[N]+'\n')
		    except KeyError:
                        R1.write("@"+thing[N].replace(" ","_")+'/1\n'+thing[S]+'\n'+'+'+'\n'+thing[A]+'\n')
			k+=1
		if " 2:" in thing[N]:
                    try:
                        match2 = db[thing[N].replace(" 2:"," 1:")]  
                    except KeyError:
                        R2.write("@"+thing[N].replace(" ","_")+'/2\n'+thing[S]+'\n'+'+'+'\n'+thing[A]+'\n')
			l+=1
                if n%1000 == 0:
                    print("Shuffleing sequence ",n)

if k == 0:
    os.remove(Singleton1)
if l == 0:
    os.remove(Singleton2)

log(Mypath, ReadNum, Filein, Fileout, Singleton1, Singleton2, n, j, k, l)

if Percent !=100:
    resamp(Percent, os.path.join(Mypath,'temp'), Mypath, Fileout, db, j/2, NumResamp)	
    print('Shuffleing done with resampleing. Check shuffle_resamp.log for a summary of the run.')
    os.remove(os.path.join(Mypath,'temp'))    
else:
    print('Shuffleing done without resampleing. Check shuffle_resamp.log for a summary of the run.')
