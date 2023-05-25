import argparse
import pandas as pd
import numpy as np
from pandas import DataFrame

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Path to the folder where the outputfile will be stored.")
argParser.add_argument("--ibdmatrix", help="Path to the folder where the input IBD genotype matrix (output from step 3) is stored.")
argParser.add_argument("--siblings", help="Path to a file with the sibling pairs you want to include in your allele frequency computations. The file should be tab delimieted, have three columns and have the header fid ID1 ID2. Each line corresponds to one sibling pair. fid: Family ID, ID1: identifier for sibling 1, ID2: identifier for sibling 2.")
argParser.add_argument("--chr", help="Chromosome")
argParser.add_argument("--trim", help="Number of SNPs to trim from the beginning and end of IBD regions for each sibling pair.")

args = argParser.parse_args()
print("args=%s" % args)

chr=args.chr
chrc="chr"+str(chr)
trim=args.trim

#inputfiles
ibdfs=args.ibdmatrix+"/trim"+str(trim)+"_"+"al1s_"+str(chr)+".csv.gz"
ibdfn=args.ibdmatrix+"/trim"+str(trim)+"_"+"al1n_"+str(chr)+".csv.gz"
#outputfile
fibd=args.out+"/trim"+str(trim)+"_"+"var1_"+str(chr)+".csv.gz"

#file with pairs to include
#The file should be tab delimieted, have three columns and have the header "fid ID1 ID2". 
#Each line corresponds to one sibling pair. fid: Family ID, ID1: identifier for sibling 1, ID2: identifier for sibling 2.
pairsb=args.siblings

print ('Input files',ibdfn,ibdfs)
print ('Sibling file', pairsb)
print ('Outputfile', fibd)
print ('Chromsome ',str(chr))
print ('Number of SNPs trimmed ',str(trim))

Gib1=pd.read_csv(ibdfn,compression='gzip')
pair=pd.read_csv(pairsb,delimiter="\t")
pair['IDa']=pair['ID1'].astype(str)+"_ibd1n_"+pair['ID2'].astype(str)
pair['IDb']=pair['ID2'].astype(str)+"_ibd1n_"+pair['ID1'].astype(str)
pairs1=list(pair['IDa'])+list(pair['IDb'])
pairs1A=list(pair['IDa'])
pairs1B=list(pair['IDb'])

snp=Gib1['SNP']
ref=Gib1['REF']
alt=Gib1['ALT']

Gib2=pd.read_csv(ibdfs,compression='gzip')
pair2=pd.read_csv(pairsb,delimiter="\t")
pair2['IDa']=pair2['ID1'].astype(str)+"_ibd1s_"+pair2['ID2'].astype(str)
pair2['IDb']=pair2['ID2'].astype(str)+"_ibd1s_"+pair2['ID1'].astype(str)
pairs2=list(pair2['IDa'])+list(pair2['IDb'])
pairs2A=list(pair2['IDa'])
pairs2B=list(pair2['IDb'])

Gib1C=np.add(Gib1[pairs1A],Gib1[pairs1B])
Gib2C=np.add(Gib2[pairs2A],Gib2[pairs2B])
GibC=np.subtract(Gib2C,Gib1C)
GibC=0.5*GibC

ibsnp=GibC.count(axis=1)
ibdiffvar=GibC.var(axis=1)

Gib1['SNP']=snp
Gib1['REF']=ref
Gib1['ALT']=alt
Gib1['ibd1']=ibsnp
Gib1['ib1diffvar']=ibdiffvar

ib=Gib1.loc[:,['SNP','REF','ALT','ibd1','ib1diffvar']]

ib.to_csv(fibd,compression='gzip',index = None, header=True)
