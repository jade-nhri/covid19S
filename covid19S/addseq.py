#!/usr/bin/env python3
import os,sys

file1=os.path.abspath(sys.argv[1])
file2=os.path.abspath(sys.argv[2])
outfile=sys.argv[3]

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
    return seq,name

d1=dict()
f1=open(file1)
with f1 as fp:
    for name,seq in read_fasta(fp):
        d1[name]=seq
f1.close()

d2=dict()
f2=open(file2)
with f2 as fp:
    for name,seq in read_fasta(fp):
        d2[name]=seq
f2.close()

for i in d1.keys():
    if i in d2.keys():
        if d1[i] in d2[i]:
            print ("{0} is identical in both".format(i))
    else:
        print ("{0} was added to outfile".format(i))
        d2[i]=d1[i]

fw=open(outfile,'w')
for i in d2.keys():
    fw.write(i+'\n')
    fw.write(d2[i]+'\n')
fw.close()
