import pandas as pd
from pandas import DataFrame
import numpy as np
import argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Outputfile")
argParser.add_argument("--ibd", help="IBD input file from snipar")
argParser.add_argument("--chr", help="Chromosome")
argParser.add_argument("--bim", help="Bim file with SNPs that are supposed to be included in the IBD matrix")
argParser.add_argument("--trim", help="How many SNPs to trim from the beginning and end of IBD regions. Default=0.")

args = argParser.parse_args()
print("args=%s" % args)

foldername=args.out
pathseg=args.ibd
chr=args.chr
trim=args.trim
pathbim=args.bim

filename1=foldername+"/trim"+str(trim)+"_ibd1_"+str(chr)+".csv.gz"
#Name of the segmentfile
pathsegfile=pathseg+"/chr_"+str(chr)+".ibd.segments.gz"
pathbimfile=pathbim+"/chr"+str(chr)+".bim"
print ('IBD segment file',pathsegfile)
print ('IBD1 output file', filename1)
print ('Number of SNPs trimmed', trim)
print ('Chromsome ',str(chr))
print ('Bimfile ',pathbimfile)

###For each chromosome, make a file with all the snps and IBD status
###The files are compressed
E = pd.read_csv(pathsegfile,compression='gzip',sep='\t')
E = pd.DataFrame(E)
E['sib']=E['ID1'].astype(str)+"_"+E['ID2'].astype(str)
E['sib']=E['sib'].astype(list)
sib=np.unique(E['sib'])
E=E.loc[E['Chr'].astype(str)==str(chr)]
E=E.loc[E['IBDType'].astype(str)=="1"]

bim = pd.read_csv(pathbimfile,sep="\t",header=None)
bim.columns=["chr","SNP","b","pos","A1","A2"]
bim=bim.loc[bim['chr'].astype(str)==str(chr)]
bim=bim.sort_values(by=['pos'])
m=bim.shape[0]

#Make a matrix with M rows (M: number of SNPs) and (S+6) columns (S: number of sibling pairs)
#The sibling pair columns indicate for which SNPs the sibling pair in question shared IBD=0 (if IBD=0 then 1 else 0)
for i in range(0,len(sib)):
        S1=E.loc[E['sib']==sib[i]]
        S1=S1.sort_values(by=['start_coordinate'])
        bim[sib[i]]='0'
        if not S1.empty:
                n=S1.shape[0]
                for j in range(0,n):
                        star=min(m,bim[bim['SNP']==S1.iloc[j,6]].index[0].astype(int)+int(trim))
                        sto=max(star,bim[bim['SNP']==S1.iloc[j,7]].index[0].astype(int)-int(trim)+1)
                        if sto>star:
                                bim[sib[i]].iloc[star:sto]=np.repeat(1,sto-star)

bim.to_csv(filename1,compression='gzip',index=None,na_rep="NA")
