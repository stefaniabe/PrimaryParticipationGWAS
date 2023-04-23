import argparse
import pandas as pd
import numpy as np
from pandas import DataFrame

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Outputfile")
argParser.add_argument("--ibdstatus", help="Path to the folder where the IBDstatus matrices is stored")
argParser.add_argument("--ibdallele", help="Path to the folder where the IBDallele matrix is stored, needs to be from analysis with trim=0")
argParser.add_argument("--chr", help="Chromosome")
argParser.add_argument("--trim", help="How many SNPs you want to trim from the beginning and end of the IBD regions")

args = argParser.parse_args()
print("args=%s" % args)

chr=args.chr
trim=args.trim
fileibd0=args.ibdstatus+"/trim"+str(trim)+"_ibd0_"+str(chr)+".csv.gz"
fileibdal0=args.ibdallele+"/trim0_ibdal0_"+str(chr)+".csv.gz"
outibd0=args.out+"/trim"+str(trim)+"_al0_"+str(chr)+".csv.gz"

print ('IBD0 status file',fileibd0)
print ('IBD0 allele file', fileibdal0)
print ('Outputfile', outibd0)
print ('Chromsome ',str(chr))
print ('Number of SNPs trimmed ',str(trim))

#Read in information about ibd regions
E = pd.read_csv(fileibd0,compression='gzip',sep=',')

#Read in information about which allele is shared
S=pd.read_csv(fileibdal0,delimiter=',',compression='gzip')
SNP=S[['SNP']]
REF=S[['REF']]
ALT=S[['ALT']]
S=S.iloc[:,3:]

e=E.shape[1]
E=E.iloc[:,6:e]
names=E.columns

names1=list(names)
names2=list(names)

for j in range(0,names.shape[0]):
        names1[j]=names[j].split("_")[0]+"_ibd0n_"+names[j].split("_")[1]
        names2[j]=names[j].split("_")[1]+"_ibd0n_"+names[j].split("_")[0]


S=S+5

S1=S[names1]
S2=S[names2]

S1=np.multiply(S1,E)
S1[S1==0]=np.nan
S1=S1-5

S2=np.multiply(S2,E)
S2[S2==0]=np.nan
S2=S2-5

S1['SNP']=SNP
S2['SNP']=SNP

S=S1.merge(S2,left_on='SNP',right_on='SNP')

S['ALT']=ALT
S['REF']=REF

###Write out compressed files
S.to_csv(outibd0,compression='gzip',index = None, header=True)
