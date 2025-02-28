#! /usr/bin/python
import argparse
import gzip
import io
from collections import Counter
from datetime import datetime

start_time = datetime.now()

def parse_user_input():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--outdir',required=True,help='Path to output directory.')
    parser.add_argument('-p','--prefix',required=True,help='Prefix for output files.')
    parser.add_argument('-r1','--read1fastq',required=True,help='Path to read 1 fastqs separated by commas.')
    parser.add_argument('-r2','--read2fastq',required=True,help='Path to read 2 fastqs searated by commas.')
    parser.add_argument('-s','--start',required=True,help='Start sequence for barcode.')
    parser.add_argument('-e','--end',required=True,help='Stop sequence for barcode.')
    parser.add_argument('-w','--whitelist',required=True,help='Path to cell barcode whitelist.')
    return parser

def check_barcode(bc):
    W=set(['A','T'])
    S=set(['G','C'])
    if bc.find('N')==-1:
        if bc[1] in S and bc[5] in S and bc[9] in S and bc[13] in S and bc[17] in S:
            if bc[3] in W and bc[7] in W and bc[11] in W and bc[15] in W and bc[19] in W:
                go=1
            else:
                go=0
        else:
            go=0
    else:
        go=0
    return go 

parser = parse_user_input()
ui = parser.parse_args()

outdir = ui.outdir
prefix = ui.prefix
r1fastq_INFILES = ui.read1fastq.split(',')
r2fastq_INFILES = ui.read2fastq.split(',')
start=ui.start
end =ui.end
lstart=len(start)
lend=len(end)
bclen=20
cbclen=16
umilen=12

whitelist = ui.whitelist
print(whitelist)
cbcset = set([line.split()[0] for line in open(whitelist)])

bcdict={}

for r1,r2 in zip(r1fastq_INFILES,r2fastq_INFILES):
    print(r1, r2)
    i=0
    with io.BufferedReader(gzip.open(r1,'rb')) as f1, io.BufferedReader(gzip.open(r2,'rb')) as f2:
        for line1,line2 in zip(f1,f2):
            if i==0:
                i+=1
            elif i==1:
                l2 = line2.decode()
                st = l2.find(start)
                if st>-1:
                    bc = l2[st+lstart:st+lstart+bclen]
                    if len(bc)==bclen and check_barcode(bc)==1:
                        l1=line1.decode()
                        cbc=l1[0:cbclen]
                        umi=l1[cbclen:cbclen+umilen]
                        bc=bc+'_'+umi
                        if cbc not in bcdict:
                            bcdict[cbc]=[bc]
                        else:
                            bcdict[cbc].append(bc)               
                else:
                    sp=l2.find(end)
                    if sp>-1:
                        bc=l2[sp-bclen:sp]
                        if len(bc)==bclen and check_barcode(bc)==1:
                            l1=line1.decode()
                            cbc=l1[0:cbclen]       
                            umi=l1[cbclen:cbclen+umilen]
                            bc=bc+'_'+umi 
                            if cbc not in bcdict:
                                bcdict[cbc]=[bc]
                            else:
                                bcdict[cbc].append(bc)     

outfile = outdir+'/'+prefix+'.tracer.txt'
print(outfile)
with open(outfile,'w') as g:
    for cbc in bcdict:
        if cbc in cbcset:
            vec1=set(bcdict[cbc])
            vec2=[bc.split('_')[0] for bc in vec1]
            cts = Counter(vec2).most_common()
            st = cbc
            for ct in cts:
                st1=ct[0]
                st2=str(ct[1])
                st=st+'\t'+st1+'\t'+st2
            st=st+'\n'
            g.write(st)
         

end_time = datetime.now()
print(f'DURATION {end_time - start_time}')
