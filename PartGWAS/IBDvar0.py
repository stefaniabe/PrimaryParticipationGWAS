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

#inputfile
ibdf=args.ibdmatrix+"/trim"+str(trim)+"_"+"al0_"+str(chr)+".csv.gz"
#outputfile
#variance file
fibd=args.out+"/trim"+str(trim)+"_"+"var0_"+str(chr)+".csv.gz"

#file with pairs to include
#The file should be tab delimieted, have three columns and have the header "fid ID1 ID2". 
#Each line corresponds to one sibling pair. fid: Family ID, ID1: identifier for sibling 1, ID2: identifier for sibling 2.
pairs=args.siblings

print ('Input file',ibdf)
print ('Sibling file', pairs)
print ('Outputfile', fibd)
print ('Chromsome ',str(chr))
print ('Number of SNPs trimmed ',str(trim))

Gib=pd.read_csv(ibdf,compression='gzip')
pair=pd.read_csv(pairs,delimiter="\t")
pair['IDa']=pair['ID1'].astype(str)+"_ibd0n_"+pair['ID2'].astype(str)
pair['IDb']=pair['ID2'].astype(str)+"_ibd0n_"+pair['ID1'].astype(str)
pairs=list(pair['IDa'])+list(pair['IDb'])
pairs1A=list(pair['IDa'])
pairs1B=list(pair['IDb'])

snp=Gib['SNP']
ref=Gib['REF']
alt=Gib['ALT']

Gib=Gib[pairs]
ibsnp=Gib.count(axis=1)

Gib['SNP']=snp
Gib['REF']=ref
Gib['ALT']=alt
Gib['ibd0']=ibsnp*0.5

#Aggregate the alleles of siblings per pair
Gib1=np.add(Gib[pairs1A],Gib[pairs1B])
#Want to have information about the allele frequency per pair
Gib1=Gib1*(0.25)
#Compute the variance per SNP
ib0meanvar=Gib1.var(axis=1)
Gib['ib0meanvar']=ib0meanvar

ib=Gib.loc[:,['SNP','REF','ALT','ibd0','ib0meanvar']]

ib.to_csv(fibd,compression='gzip',index = None, header=True)
