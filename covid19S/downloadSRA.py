#!/usr/bin/env python3 
import os, subprocess
import sys
import pandas as pd 
import multiprocessing as mp
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-o', help='an output folder')
parser.add_argument('-t', help='threads (default=100)')
parser.add_argument('--project', help='SRA projcet ID')
parser.add_argument('--nanopore', help='only download Nanopore sequences')

args = parser.parse_args() 
t=100 
only_nanopore=False
argv=sys.argv 
if '-o' in argv: 
    outdir=argv[argv.index('-o')+1] 
if '-t' in argv: 
    t=int(argv[argv.index('-t')+1]) 
if '--project' in argv:
    sraprojectID=argv[argv.index('--project')+1]
if '--nanopore' in argv:
    only_nanopore=True


def fqdump(RunID): 
    comm='fastq-dump {0}'.format(RunID) 
    print (comm) 
    subprocess.getoutput(comm)

cwd=os.getcwd() 
outdir=argv[argv.index('-o')+1] 
if  not os.path.exists(outdir): 
    os.mkdir(outdir) 
os.chdir(outdir) 


comm="wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term={0}' -O - > {1}".format(sraprojectID,'sralist.csv') 
print (comm) 
subprocess.getoutput(comm) 
df=pd.read_csv('sralist.csv')
platforms=np.unique(df['Platform'].values)
print (platforms)

platform_name=''
for i in platforms:
    tmp=i.lower()
    #print (tmp)
    if 'nanopore' in tmp:
        platform_name=i
        print(platform_name)

print (df) 
if only_nanopore and platform_name!='':
    df=df[df['Platform']==platform_name]
print (df)
runlist=df['Run'].values 
print (runlist) 
def main(): 
    po=mp.Pool(t) 
    for i in runlist: 
        po.apply_async(fqdump,args=[i]) 
    po.close() 
    po.join() 
if __name__=='__main__': 
    main()


