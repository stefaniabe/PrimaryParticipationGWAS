import argparse
import pandas as pd
import numpy as np
from pandas import DataFrame

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Outputfolder")
argParser.add_argument("--ibdstatus", help="Path to the folder where the IBDstatus matrix is stored, needs to be from analysis with trim=0")
argParser.add_argument("--chr", help="Chromosome")
argParser.add_argument("--bim", help="Path to a bimfile with the SNPs in the vcf files.")
argParser.add_argument("--vcf", help="Path to the chromosome specific vcf files with the phased genotypes of individuals with a sibling in the dataset (compressed file ending with .gz)")
argParser.add_argument("--headerline", help="In which line is the header in the vcf files.")

args = argParser.parse_args()
print("args=%s" % args)

chr=args.chr
#Run this with trim=0. Only need to run this script once
trim=0

#vcfheaderline

vcfheader=int(args.headerline)-1

#Chromosome specific vcf file
#The vcf file has the header in row nr. 4
pathbed=args.vcf+"/chr"+str(chr)+".vcf.gz"

foldername=args.out
pathibd=args.ibdstatus+"/trim"+str(trim)+"_ibd0_"+str(chr)+".csv.gz"

#output file
fileibd0=foldername+"/trim"+str(trim)+"_ibdal0_"+str(chr)+".csv.gz"

#bim file
pathbim=args.bim
pathbimfile=pathbim+"/chr"+str(chr)+".bim"

print ('IBD status file',pathibd)
print ('IBD0 output file', fileibd0)
print ('VCF file', pathbed)
print ('Chromsome ',str(chr))
print ('Bimfile ',pathbimfile)

#Read in information about ibd regions to make sure we have the siblings in the same order
#Just need the header
E = pd.read_csv(pathibd,compression='gzip',sep=',',nrows=0)
sib=E.columns
sib=sib[6:sib.shape[0]]

G=pd.read_csv(pathbed,delimiter='\t',skiprows=vcfheader,compression='gzip')
#Make sure that the columns SNP name and CHR are named that in the dataframe
G.rename(columns={'ID': 'SNP', '#CHROM': 'CHR'}, inplace=True)
#Quick fix in case the SNP column is on the format CHR,SNP.
ccc=","+str(chr)
G['SNP']=G['SNP'].str.replace(ccc,"")

bim = pd.read_csv(pathbimfile,sep="\t",header=None)
bim.columns=["chr","SNP","b","pos","A1","A2"]
bim=bim.loc[bim['chr'].astype(str)==str(chr)]

G=G.merge(bim[['SNP']],left_on='SNP',right_on='SNP',how='inner')
G=G.sort_values(by=['POS'])

Eibd0=G[['SNP','REF','ALT']]

for i in range(0,len(sib)):
        s=sib[i]
        s1=s[:7]
        s2=s[8:]
        S=G[['SNP','POS',s1,s2]]
        s1ibd0=s1+"_ibd0n_"+s2
        s2ibd0=s2+"_ibd0n_"+s1
        S[s1ibd0]=(S[s1].str[0].astype(int)+S[s1].str[2].astype(int))
        S[s2ibd0]=(S[s2].str[0].astype(int)+S[s2].str[2].astype(int))
        Eibd0=Eibd0.merge(S[['SNP',s1ibd0,s2ibd0]],left_on='SNP',right_on='SNP',how='outer')

###Write out a compressed file
Eibd0.to_csv(fileibd0,compression='gzip',index = None, header=True)
