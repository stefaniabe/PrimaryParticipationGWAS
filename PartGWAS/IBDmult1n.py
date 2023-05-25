import argparse
import pandas as pd
import numpy as np
from pandas import DataFrame

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Outputfile")
argParser.add_argument("--ibdstatus", help="Path to the folder where the IBDstatus matrices is stored")
argParser.add_argument("--ibdallele", help="Path to the folder where the IBDallele matrix is stored, needs to be from analysis with trim=0")
argParser.add_argument("--chr", help="Chromosome")
argParser.add_argument("--maf", help="Path to the folder where the maf information is stored. Note, all the SNPs in the IBDstatus matrix and the IBDallele matrix need to be in the maf file.")
argParser.add_argument("--trim", help="How many SNPs you want to trim from the beginning and end of the IBD regions")

args = argParser.parse_args()
print("args=%s" % args)

foldername=args.out
chr=args.chr
trim=args.trim
fileibd1=args.ibdstatus+"/trim"+str(trim)+"_ibd1_"+str(chr)+".csv.gz"
fileibdal1=args.ibdallele+"/trim0_ibdal1n_"+str(chr)+".csv.gz"
mafpath=args.maf+"/chr"+str(chr)+".frq"
outibd1=args.out+"/trim"+str(trim)+"_al1n_"+str(chr)+".csv.gz"

print ('IBD1 status file',fileibd1)
print ('IBD1 not-shared allele file', fileibdal1)
print ('Outputfile', outibd1)
print ('Chromsome ',str(chr))
print ('MAF file ',mafpath)
print ('Number of SNPs trimmed ',str(trim))

#Read in information about ibd regions
E = pd.read_csv(fileibd1,compression='gzip',sep=',')

#Read in information about which allele is shared
S=pd.read_csv(fileibdal1,delimiter=',',compression='gzip')
SNP=S[['SNP']]
REF=S[['REF']]
ALT=S[['ALT']]
#Read in maf information
#Use maf information to impute the "shared" allele when we are unable to infer the shared allele from phasing information
maf=pd.read_csv(mafpath,delimiter="\t")
maf=maf.merge(S[['SNP','REF']],left_on='SNP',right_on='SNP')
maf['meang']=np.where(maf['A1']==maf['REF'],1-maf['MAF'],maf['MAF'])

S=S.merge(maf[['SNP','meang']],left_on='SNP',right_on='SNP')
s=S.shape[1]-1
meang=S[['meang']]
S=S.iloc[:,3:s]

for j in range(0,S.shape[0]):
        S.iloc[j,:]=np.where(abs(S.iloc[j,:])>500,meang.iloc[j],S.iloc[j,:])

e=E.shape[1]
E=E.iloc[:,6:e]
names=E.columns

names1=list(names)
names2=list(names)

for j in range(0,names.shape[0]):
        names1[j]=names[j].split("_")[0]+"_ibd1n_"+names[j].split("_")[1]
        names2[j]=names[j].split("_")[1]+"_ibd1n_"+names[j].split("_")[0]
        
        
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
S.to_csv(outibd1,compression='gzip',index = None, header=True)
