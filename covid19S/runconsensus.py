#!/usr/bin/env python3
import os,sys
import subprocess
import argparse
import pysam
import csv


parser = argparse.ArgumentParser()
parser.add_argument('-i', help='an input folder')
parser.add_argument('-o', help='an output folder')
parser.add_argument('-t', help='threads (default=100)')
parser.add_argument('-s', help='specify a sample name')
parser.add_argument('-m', help='min length of consensus (default=3000 bp)')
parser.add_argument('-r', help='the path to the reference sequence')


args = parser.parse_args()

minlength=3000
t=100
bc=''
ref='/MyData/mnt/nas40T/jade/202102/SRA/Sgene.fasta'

argv=sys.argv
if '-i' in argv:
    indir=argv[argv.index('-i')+1]
if '-o' in argv:
    outdir=argv[argv.index('-o')+1]
if '-t' in argv: 
    t=int(argv[argv.index('-t')+1]) 
if '-s' in argv:
    bc=argv[argv.index('-s')+1]
if '-m' in argv:
    minlength=int(argv[argv.index('-m')+1])
if '-r' in argv:
    ref=argv[argv.index('-r')+1]

indir=os.path.abspath(indir)
outdir=os.path.abspath(outdir)
ref=os.path.abspath(ref)

myfiles=[x for x in os.listdir(indir) if '.fastq' in x]
#print (myfiles)

if not os.path.exists(outdir):
    os.mkdir(outdir)

def read_fasta(fp): 
    name, seq = None, [] 
    for line in fp: 
        line = line.rstrip() 
        if line.startswith(">"): 
            if name: yield (name, ''.join(seq)) 
            name, seq = line, [] 
        else: 
            seq.append(line) 
    if name: yield (name, ''.join(seq)) 


def getdepth(indir):
    indir=os.path.abspath(indir)
    tmppath=os.path.split(indir)
    infa=ref
    inbam=os.path.join(indir,'calls_to_draft.bam')
    if os.path.exists(inbam) and os.path.getsize(inbam)>0:
        bamfile=pysam.AlignmentFile(inbam,'rb')
    else:
        return
    d=dict()
    f=open(infa)
    with f as fp:
        for name, seq in read_fasta(fp):
            d[name]=seq
    f.close()
    dout=dict()
    allpos=set()
    for i in range(1,len(seq)+1):
        allpos.add(i)
    for i in d.keys():
        #print (i)
        pos=[]
        posset=set()
        for pileupcolumn in bamfile.pileup(bamfile.get_reference_name(0),0,):
            missingdepth=100-pileupcolumn.n
            posset.add(pileupcolumn.pos+1)
            if missingdepth>0:
                #print ('{0}\t{1}'.format(pileupcolumn.pos,missingdepth))
                dout[pileupcolumn.pos+1]=missingdepth
        zeropos=allpos-posset
        for j in zeropos:
            dout[j]=100
    with open(os.path.join(tmppath[0],tmppath[1]+'_md.txt'),'w') as f:
        for key in dout.keys():
            f.write('{0}\t{1}\n'.format(key,dout[key]))
    return (dout)

myfile=[]
for fqfile in myfiles:
    if bc!='' and bc in fqfile:
        print (fqfile)
        myfile.append(fqfile)
    if bc=='':
        myfile.append(fqfile)
#print (myfile)

fw=open(os.path.join(outdir,'consensus.fasta'),'w')
for fqfile in sorted(myfile):
    out=fqfile.replace('.fastq','')
    #print (out)
    if os.path.exists(ref) and os.path.getsize(ref)>0:
        if os.path.exists(os.path.join(outdir,out+'_final_old.fa')):
            continue
        if os.path.exists(os.path.join(outdir,out+'_final.fa')):
            dout=dict()
            ff=open(os.path.join(outdir,out+'_final.fa'))
            nseq=0
            with ff as fp:
                for namef,seqf in read_fasta(fp):
                    nseq+=1
                    dout[out]=seqf
            ff.close()
            tempfiles=[x for x in os.listdir(outdir) if '_ref_5.fa' in x]
            #print (tempfiles)
            #print (len(dout.keys()))
            if nseq==1 and out+'_ref_5.fa' not in tempfiles and len(seqf)>minlength:
                fw.write('>'+out+'\n')
                fw.write(dout[out]+'\n')
                #print (out)

            else:
                print('{0} was not included in the consensus.fasta'.format(out))
                comm='mv {0} {1}'.format(os.path.join(outdir,out+'_final.fa'),os.path.join(outdir,out+'_final.fa').replace('_final.fa','_final_old.fa'))
                print (comm)
                subprocess.getoutput(comm)


            continue
        comm='medaka_consensus -i {2}/{0} -d {3} -o {4}/{1} -g -t 100 -f'.format(fqfile,out,indir,ref,outdir)
        #print (comm)
        stdout=subprocess.getoutput(comm)
        plotd=getdepth(os.path.join(outdir,out))
        #print (plotd)
        for i in range(1,6):
            if not os.path.exists('{0}/{1}/consensus.fasta'.format(outdir,out)) or os.path.getsize('{0}/{1}/consensus.fasta'.format(outdir,out))==0:
                print ('No consensus sequence')
                break
            comm='cp {2}/{0}/consensus.fasta {2}/{0}_ref_{1}.fa'.format(out,i,outdir)
            #print (comm)
            subprocess.getoutput(comm)
            #comm='rm {0} -rf'.format(out)
            #print (comm)
            #subprocess.getoutput(comm)
            #print (out+'_ref_'+str(i)+'.fa')
            if not os.path.exists(os.path.join(outdir,out+'_ref_'+str(i)+'.fa')) or os.path.getsize(os.path.join(outdir,out+'_ref_'+str(i)+'.fa'))==0:
                break
            comm='medaka_consensus -i {3}/{0} -d {4}/{1}_ref_{2}.fa -o {4}/{1} -g -t 100 -f'.format(fqfile,out,i,indir,outdir)
            #print (comm)
            stdout=subprocess.getoutput(comm)
            #print (stdout)
            f1=open(os.path.join(outdir,out+'_ref_'+str(i)+'.fa'))
            d1=dict()
            with f1 as fp:
                for name1, seq1 in read_fasta(fp):
                    #print (seq1)
                    sseq1=''.join(seq1)
            #print (sseq1)
            f1.close()
            f2=open(os.path.join(outdir,out+'/consensus.fasta'))
            d2=dict()
            with f2 as fp:
                for name2, seq2 in read_fasta(fp):
                    #print (seq2)
                    sseq2=''.join(seq2)
            f2.close()

            if sseq1==sseq2:
                break

    if os.path.exists('{0}/{1}_ref_{2}.fa'.format(outdir,out,i)) and os.path.getsize('{0}/{1}_ref_{2}.fa'.format(outdir,out,i))>0:
        comm='cp {2}/{0}_ref_{1}.fa {2}/{0}_final.fa'.format(out,i,outdir)
        #print (comm)
        subprocess.getoutput(comm)
        ff=open(os.path.join(outdir,out+'_final.fa'))
        dout=dict()
        nseq=0
        with ff as fp:
            for namef,seqf in read_fasta(fp):
                nseq+=1
                dout[out]=seqf
        ff.close()
        tempfiles=[x for x in os.listdir(outdir) if '_ref_5.fa' in x]
        #print (tempfiles)
        if nseq==1 and out+'_ref_5.fa' not in tempfiles:
            fw.write('>'+out+'\n')
            fw.write(dout[out]+'\n')

        else:
             print('{0} was not included in the consensus.fasta'.format(out))
             comm='mv {0} {1}'.format(os.path.join(outdir,out+'_final.fa'),os.path.join(outdir,out+'_final.fa').replace('_final.fa','_final_old.fa'))
             #print (comm)
             subprocess.getoutput(comm)



fw.close()





