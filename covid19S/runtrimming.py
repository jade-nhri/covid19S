#!/usr/bin/env python3 
import os,sys 
import subprocess 
import argparse
import multiprocessing as mp

parser = argparse.ArgumentParser() 
parser.add_argument('-i', help='an input folder') 
parser.add_argument('-o', help='an output folder') 
parser.add_argument('-t', help='threads (default=100)')
parser.add_argument('-l', help='trimming length (default=100 bp)')
args = parser.parse_args() 
t=100
tlen=100

argv=sys.argv
if '-i' in argv:
    indir=argv[argv.index('-i')+1]
if '-o' in argv:
    outdir=argv[argv.index('-o')+1]
if '-t' in argv:
    t=int(argv[argv.index('-t')+1])
if '-l' in argv:
    tlen=int(argv[argv.index('-l')+1])

indir=os.path.abspath(indir) 
outdir=os.path.abspath(outdir)

def trimfq(infile):
    comm='cat {0}/{1} | seqkit subseq -r {2}:-{2} > {3}/{1}'.format(indir,infile,tlen,outdir)
    print (comm)
    subprocess.getoutput(comm)

os.chdir(indir)
myfiles=[x for x in os.listdir() if '.fastq' in x]
if not os.path.exists(outdir):
    os.mkdir(outdir)
def main():
    po=mp.Pool(t)
    for i in myfiles:
        po.apply_async(trimfq,args=[i])
    po.close()
    po.join()

if __name__=='__main__':
    main()
