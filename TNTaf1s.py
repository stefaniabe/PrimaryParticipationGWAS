import argparse
import pandas as pd
import numpy as np
from pandas import DataFrame

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Path to the folder where the outputfile will be stored.")
argParser.add_argument("--ibdmatrix", help="Path to the folder where the input IBD genotype matrix (output from step 1) is stored.")
argParser.add_argument("--po_pairs", help="Path to a file with the parent-offspring pairs you want to include in your allele frequency computations. The file should be tab delimieted, have three columns and have the header fid ID1 ID2. Each line corresponds to one parent-offspring pair. fid: Family ID, ID1: identifier for the offspring, ID2: identifier for the parent.")
argParser.add_argument("--chr", help="Chromosome")

args = argParser.parse_args()
print("args=%s" % args)

chr=args.chr
chrc="chr"+str(chr)

#inputfile
ibdf=args.ibdmatrix+"/TNT_al1s_"+str(chr)+".csv.gz"
#outputfile
#frequency file
fibd=args.out+"/TNT_af1s_"+str(chr)+".csv.gz"

#file with pairs to include
#The file should be tab delimieted, have three columns and have the header "fid ID1 ID2". 
#Each line corresponds to one parent-offspring pair. fid: Family ID, ID1: identifier for sibling 1, ID2: identifier for sibling 2.
pairs=args.po_pairs

print ('Chromosome ',str(chr))
print ('Input file',ibdf)
print ('Parent-offspring file', pairs)
print ('Outputfile', fibd)

Gib=pd.read_csv(ibdf,compression='gzip')
pair=pd.read_csv(pairs,delimiter="\t")
#I just want to examine the non-transmitted vs transmitted alleles of parents
pair['IDa']=pair['ID2'].astype(str)+"_ibd1s_"+pair['ID1'].astype(str)
pairs=list(pair['IDa'])

snp=Gib['SNP']
ref=Gib['REF']
alt=Gib['ALT']

Gib=Gib[pairs]

ibsnp=Gib.count(axis=1)
ibmean=Gib.mean(axis=1)

Gib['SNP']=snp
Gib['REF']=ref
Gib['ALT']=alt
Gib['ibd1s']=ibsnp
Gib['ib1smean']=ibmean

ib=Gib.loc[:,['SNP','REF','ALT','ibd1s','ib1smean']]

ib.to_csv(fibd,compression='gzip',index = None, header=True)
