#!/usr/bin/env python3 
import os,sys
import subprocess
import edlib 
from Bio.Seq import Seq 
import pysam 
from collections import Counter 
import argparse
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='an input folder')
parser.add_argument('-o', help='an output folder')
parser.add_argument('-r', help='the path to the reference sequence')
parser.add_argument('-q', help='an input folder containing reads in fastq')


refpath='/opt/covid19S/covid19S/Sgene.fasta'
refpathp='/opt/covid19S/covid19S/Sprotein.fa'


argv=sys.argv
if '-i' in argv:
    indir=argv[argv.index('-i')+1]
if '-o' in argv:
    outdir=argv[argv.index('-o')+1]
if '-r' in argv:
    refpath=argv[argv.index('-r')+1]
if '-q' in argv:
    fqdir=argv[argv.index('-q')+1]
    fqdir=os.path.abspath(fqdir)
if '-rp' in argv:
    refpathp=argv[argv.index('-rp')+1]

indir=os.path.abspath(indir)
outdir=os.path.abspath(outdir)

if not os.path.exists(outdir):
    os.mkdir(outdir)


outfile=os.path.join(outdir,'output.txt')
outfa=os.path.join(outdir,'output.fasta')

os.chdir(indir)
#myfiles=[x for x in os.listdir() if '_final' in x and ('ERR4365241_final.fa' in x or 'ERR4867938_final.fa' in x or 'ERR4890482_final.fa' in x)] 
myfiles=[x for x in os.listdir() if '_final.fa' in x]

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
def processgap (query,subject): 
    gappos=0 
    k=4
    if '-' in query: 
        for i in range(len(query)): 
            if '-AA' in query[i:i+3] or 'AA-' in query[i:i+3]: 
                gappos=i 
                if subject[gappos:gappos+3]=='AAA': 
                    print ("  gap in query:"+query[gappos-k:gappos+3+k]+"==>"+query[gappos-k:gappos]+'AAA'+query[gappos+3:gappos+3+k]) 
                    query=query.replace(query[gappos-k:gappos+3+k],query[gappos-k:gappos]+'AAA'+query[gappos+3:gappos+3+k]) 
            if '-TT' in query[i:i+3] or 'TT-' in query[i:i+3]: 
                gappos=i 
                if subject[gappos:gappos+3]=='TTT': 
                    print ("  gap in query:"+query[gappos-k:gappos+3+k]+"==>"+query[gappos-k:gappos]+'TTT'+query[gappos+3:gappos+3+k]) 
                    query=query.replace(query[gappos-k:gappos+3+k],query[gappos-k:gappos]+'TTT'+query[gappos+3:gappos+3+k]) 
            if '-CC' in query[i:i+3] or 'CC-' in query[i:i+3]: 
                gappos=i 
                if subject[gappos:gappos+3]=='CCC': 
                    print ("  gap in query:"+query[gappos-k:gappos+3+k]+"==>"+query[gappos-k:gappos]+'CCC'+query[gappos+3:gappos+3+k]) 
                    query=query.replace(query[gappos-k:gappos+3+k],query[gappos-k:gappos]+'CCC'+query[gappos+3:gappos+3+k]) 
            if '-GG' in query[i:i+3] or 'GG-' in query[i:i+3]: 
                gappos=i 
                if subject[gappos:gappos+3]=='GGG': 
                    print ("  gap in query:"+query[gappos-k:gappos+3+k]+"==>"+query[gappos-k:gappos]+'GGG'+query[gappos+3:gappos+3+k]) 
                    query=query.replace(query[gappos-k:gappos+3+k],query[gappos-k:gappos]+'GGG'+query[gappos+3:gappos+3+k]) 
    if '-' in subject: 
        for i in range(len(subject)): 
            if '-AA' in subject[i:i+3] or 'AA-' in subject[i:i+3]: 
                gappos=i 
                if query[gappos:gappos+3]=='AAA': 
                    print ("  gap in subject:"+query[gappos-k:gappos+3+k]+"==>"+query[gappos-k:gappos]+'AA'+query[gappos+3:gappos+3+k]) 
                    query=query.replace(query[gappos-k:gappos+3+k],query[gappos-k:gappos]+'AA'+query[gappos+3:gappos+3+k]) 
                    subject=subject.replace(subject[gappos-k:gappos+3+k],subject[gappos-k:gappos+3+k].replace('-','')) 
            if '-TT' in subject[i:i+3] or 'TT-' in subject[i:i+3]: 
                gappos=i 
                if query[gappos:gappos+3]=='TTT': 
                    print ("  gap in subject:"+query[gappos-k:gappos+3+k]+"==>"+query[gappos-k:gappos]+'TT'+query[gappos+3:gappos+3+k]) 
                    query=query.replace(query[gappos-k:gappos+3+k],query[gappos-k:gappos]+'TT'+query[gappos+3:gappos+3+k]) 
                    subject=subject.replace(subject[gappos-k:gappos+3+k],subject[gappos-k:gappos+3+k].replace('-','')) 
            if '-CC' in subject[i:i+3] or 'CC-' in subject[i:i+3]: 
                gappos=i 
                if query[gappos:gappos+3]=='CCC': 
                    print ("  gap in subject:"+query[gappos-k:gappos+3+k]+"==>"+query[gappos-k:gappos]+'CC'+query[gappos+3:gappos+3+k]) 
                    query=query.replace(query[gappos-k:gappos+3+k],query[gappos-k:gappos]+'CC'+query[gappos+3:gappos+3+k]) 
                    subject=subject.replace(subject[gappos-k:gappos+3+k],subject[gappos-k:gappos+3+k].replace('-','')) 
            if '-GG' in subject[i:i+3] or 'GG-' in subject[i:i+3]: 
                gappos=i 
                if query[gappos:gappos+3]=='GGG': 
                    print ("  gap in subject:"+query[gappos-k:gappos+3+k]+"==>"+query[gappos-k:gappos]+'GG'+query[gappos+3:gappos+3+k]) 
                    query=query.replace(query[gappos-k:gappos+3+k],query[gappos-k:gappos]+'GG'+query[gappos+3:gappos+3+k]) 
                    subject=subject.replace(subject[gappos-k:gappos+3+k],subject[gappos-k:gappos+3+k].replace('-','')) 
    return (query) 

def countstar(seq):
    N=list(seq).count('*')
    return (N)
def countgap(seq):
    N=list(seq).count('-')
    return (N)


def trim4trans(inseq,normal):
    try:
        coding_dna=Seq(inseq)
        cds=coding_dna.translate()
        if countstar(cds)<=1:
            #print ('0')
            normal=True
            return (inseq,normal)
        else:
            coding_dna=Seq(inseq[1:])
            cds=coding_dna.translate()
            if countstar(cds)<=1:
                #print ('1')
                normal=True
                return(inseq[1:],normal)
            else:
                coding_dna=Seq(inseq[2:])
                cds=coding_dna.translate()
                if countstar(cds)<=1:
                    #print ('2')
                    normal=True
                    return(inseq[2:],normal)
                else:
                    #print ('False')
                    return(inseq,False)
    except:
        #print ('Oops')
        return (inseq,False)


    

def getconf(indir):
    indir=os.path.abspath(indir)
    infa=os.path.join(indir,'consensus.fasta')
    inbam=os.path.join(indir,'calls_to_draft.bam')
    bamfile=pysam.AlignmentFile(inbam,'rb')
    d=dict()
    f=open(infa)
    with f as fp:
        for name, seq in read_fasta(fp):
            d[name]=seq
    f.close()
    dout=dict()
    seqout=list(seq)
    for i in d.keys():
        #print (i)
        pos=[]
        for pileupcolumn in bamfile.pileup(bamfile.get_reference_name(0),0,):
            #print ('{0}\t{1}'.format(pileupcolumn.pos,pileupcolumn.n))
            if pileupcolumn.n>=50:
                pos.append(pileupcolumn.pos-1)
    for j in range(0,len(seqout)):
        if j not in pos:
            seqout[j]='N'
    tmpseq=''.join(seqout)
    #print (tmpseq)
    temp=tmpseq.split('N')
    #print (temp)
    nseq=0
    for seq in temp:
        if seq!='':
            nseq+=1
            #print (seq)
            dout['Seq'+str(nseq)]=seq
    #print (dout)
    return(dout)

fref=open(refpathp)
with fref as fp:
    for name,seq in read_fasta(fp):
        target=seq
fref.close()
#print (target)
fref=open(refpath)
with fref as fp:
    for name,seq in read_fasta(fp):
        subject=seq
fref.close()

mylist=[]
redofiles=[]
fw=open(outfile,'w')
faw=open(outfa,'w')
fwcsv=open(os.path.join(outdir,'Result.csv'),'w') 
fwcsv.write('Name,Variants,NT,AA\n')
for i in sorted(myfiles):
    comm='grep ">" {0} | wc -l'.format(i)
    #print (comm)
    Nseq=subprocess.getoutput(comm)
    
    d=dict()
    filename=i.split('_')[0]
    mylist.append(filename)
    if int(Nseq)>1:
        print (i)
        fw.write('{0}\t{1}\n'.format(filename,'multiple-fragment consensus'))
        continue
    f=open(i)
    with f as fp:
        for name, seq in read_fasta(fp):
            (seqout,trans)=trim4trans(seq,False)
            if trans:
                d[name]=seqout
            if not trans:
                alignment=edlib.align(seq,subject,task='path')
                alignment_nice=edlib.getNiceAlignment(alignment,seq,subject) 
                #print (filename)
                #print (alignment_nice['query_aligned'])
                #print (alignment_nice['target_aligned'])
                if countgap(alignment_nice['query_aligned'])<10 and countgap(alignment_nice['target_aligned'])<10:
                    corseq=processgap(alignment_nice['query_aligned'],alignment_nice['target_aligned'])
                    (seqout,trans)=trim4trans(corseq,False)
                    if trans:
                        print ('{0} was corrected'.format(i))
                        d[name]=seqout
                    else:
                        print ('{0} can not be corrected, try to checck sequencing depth...'.format(i))
                        if os.path.exists(filename):
                            d=getconf(filename)
                            for j in d.keys():
                                (seqout,trans)=trim4trans(d[j],False)
                                d[j]=seqout
                else:
                    #print ('{0} can not be corrected'.format(i))
                    if os.path.exists(filename):
                        d=getconf(filename)
                        for j in d.keys():
                            (seqout,trans)=trim4trans(d[j],False)
                            d[j]=seqout

            
    f.close()
    if d=={}:
        fw.write('{0}\t{1}\n'.format(filename,'Segment-missing amplicon'))
        fwcsv.write('{0},Segment-missing amplicon,,\n'.format(filename))

    if len(d.keys())>1:
        fw.write('{0}\t{1}\n'.format(filename,'Segment-missing amplicon'))
        fwcsv.write('{0},Segment-missing amplicon,,\n'.format(filename))
        #print ('{0}\t{1}'.format(filename,'Segment-missing amplicon'))
    if len(d.keys())==1 and trans:
        for name in d.keys():
            seqout=d[name]
        coding_dna=Seq(seqout)
        cds=coding_dna.translate()
        result=edlib.align(cds,target,task='path')
        nice=edlib.getNiceAlignment(result,cds,target)
        cdsout=nice['query_aligned']
        #print (cdsout)
        site=0
        var=[]
        for j,k in zip(target,cdsout):
            site+=1
            if j!=k and k!='X':
                #print ('{0}{1}{2}'.format(j,site,k))
                var.append('{0}{1}{2}'.format(j,site,k))
        if var==[]:
            varout=''
        else:
            varout=';'.join(var)
        Nvar=len(var)
        #print ('{0}\t{1}'.format(filename,varout))
        if Nvar<20:
            fw.write('{0}\t{1}\n'.format(filename,varout))
            faw.write('>{0}\n'.format(filename))
            faw.write('{0}\n'.format(seqout))
            fwcsv.write('{0},{1},{2},{3}\n'.format(filename,varout,seqout,cdsout)) 
        else:
            redofiles.append(i)

    if len(d.keys())==1 and not trans:
        print ('{0}\tsingle but not translated'.format(filename))
        fw.write('{0}\t{1}\n'.format(filename,'Segment-missing amplicon'))
        fwcsv.write('{0},multiple-fragment consensus,,\n'.format(filename))

#print (redofiles)

for i in redofiles:
    d=dict()
    filename=i.split('_')[0]
    f=open(i)
    with f as fp:
        for name, seq in read_fasta(fp):
            d[name]=seq
    f.close()
    d=getconf(filename)
    if len(d.keys())>1:
        fw.write('{0}\t{1}\n'.format(filename,'Segment-missing amplicon'))
        fwcsv.write('{0},Segment-missing amplicon,,\n'.format(filename))

    else:
        seq=d.values()
        if not d or len(seq)<2000:
            print ('{0} is too short'.format(filename))
            fw.write('{0}\t{1}\n'.format(filename,'Segment-missing amplicon'))
            fwcsv.write('{0},Segment-missing amplicon,,\n'.format(filename))
        else:
            for l in d.keys():
                seq=d[l]
                (seqout,trans)=trim4trans(seq,False)
                if trans:
                    print ('Check {0}'.format(filename))
                    d[l]=seqout
                    coding_dna=Seq(seqout)
                    cds=coding_dna.translate()
                    result=edlib.align(cds,target,task='path')
                    nice=edlib.getNiceAlignment(result,cds,target)
                    cdsout=nice['query_aligned']
                    site=0
                    var=[]
                    for j,k in zip(target,cdsout):
                        site+=1
                        if j!=k and k!='X':
                            #print ('{0}{1}{2}'.format(j,site,k))
                            var.append('{0}{1}{2}'.format(j,site,k))
                    if var==[]:
                        varout=''
                    else:
                        varout=';'.join(var)
                fw.write('{0}\t{1}\n'.format(filename,varout))
                faw.write('>{0}\n'.format(filename))
                faw.write('{0}\n'.format(seqout))
                fwcsv.write('{0},{1},{2},{3}\n'.format(filename,varout,seqout,cdsout))


if '-q' in argv:
    #print (fqdir)
    os.chdir(fqdir)
    myfqs=[x for x in os.listdir() if '.fastq' in x]
    for i in myfqs:
        filename=i.replace('.fastq','')
        if filename not in mylist:
            fw.write('{0}\tmultiple-fragment consensus\n'.format(filename))
            fwcsv.write('{0},multiple-fragment consensus,,\n'.format(filename))

fwcsv.close()
faw.close()
fw.close()


