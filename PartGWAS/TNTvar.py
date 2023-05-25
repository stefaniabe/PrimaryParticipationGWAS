import argparse
import pandas as pd
import numpy as np
from pandas import DataFrame

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Path to the folder where the outputfile will be stored.")
argParser.add_argument("--ibdmatrix", help="Path to the folder where the input IBD genotype matrix (output from the previous step) is stored.")
argParser.add_argument("--po_pairs", help="Path to a file with the parent-offspring pairs you want to include in your allele frequency computations. The file should be tab delimieted, have three columns and have the header fid ID1 ID2. Each line corresponds to one parent-offspring pair. fid: Family ID, ID1: identifier for the offspring, ID2: identifier for the parent.")
argParser.add_argument("--chr", help="Chromosome")

args = argParser.parse_args()
print("args=%s" % args)

chr=args.chr
chrc="chr"+str(chr)

#inputfiles
ibdfn=args.ibdmatrix+"/TNT_al1n_"+str(chr)+".csv.gz"
ibdfs=args.ibdmatrix+"/TNT_al1s_"+str(chr)+".csv.gz"
#outputfile
#variance file
fibd=args.out+"/TNT_var_"+str(chr)+".csv.gz"

#file with pairs to include
#The file should be tab delimieted, have three columns and have the header "fid ID1 ID2". 
#Each line corresponds to one parent-offspring pair. fid: Family ID, ID1: identifier for sibling 1, ID2: identifier for sibling 2.
pairsb=args.po_pairs

print ('Chromosome ',str(chr))
print ('Input files',ibdfs,ibdfn)
print ('Parent-offspring file', pairsb)
print ('Outputfile', fibd)

Gib1=pd.read_csv(ibdfn,compression='gzip')
Gib2=pd.read_csv(ibdfs,compression='gzip')

pair=pd.read_csv(pairsb,delimiter="\t")
pair['IDb']=pair['ID2'].astype(str)+"_ibd1n_"+pair['ID1'].astype(str)
pairs1B=list(pair['IDb'])
pair['IDa']=pair['ID2'].astype(str)+"_ibd1s_"+pair['ID1'].astype(str)
pairs2B=list(pair['IDa'])

snp=Gib1['SNP']
ref=Gib1['REF']
alt=Gib1['ALT']

Gib1C=Gib1[pairs1B]
Gib2C=Gib2[pairs2B]
GibC=np.subtract(Gib2C,Gib1C)

ibsnp=GibC.count(axis=1)
ibdiffvar=GibC.var(axis=1)

Gib1['SNP']=snp
Gib1['REF']=ref
Gib1['ALT']=alt
Gib1['ibd1']=ibsnp
Gib1['ib1diffvar']=ibdiffvar

ib=Gib1.loc[:,['SNP','REF','ALT','ibd1','ib1diffvar']]

ib.to_csv(fibd,compression='gzip',index = None, header=True)
