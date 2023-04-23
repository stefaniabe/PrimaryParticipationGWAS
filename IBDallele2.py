import argparse
import pandas as pd
import numpy as np
from pandas import DataFrame

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Outputfile")
argParser.add_argument("--ibdstatus", help="Path to the folder where the IBDstatus matrix is stored, needs to be from analysis with trim=0")
argParser.add_argument("--chr", help="Chromosome")
argParser.add_argument("--bim", help="Path to a bimfile with the SNPs in the vcf files that you want to include in the IBD matrices.")
argParser.add_argument("--vcf", help="Path to the chromosome specific vcf files with the phased genotypes of individuals with a sibling in the dataset (compressed file ending with .gz)")
argParser.add_argument("--headerline", help="In which line is the header in the vcf files.")

args = argParser.parse_args()
print("args=%s" % args)

foldername=args.out
chr=args.chr

vcfheader=int(args.headerline)-1

#Run this with trim=0. Only need to run this script once
trim=0

#File with chromosome specific vcf files
#The vcf file has the header in row nr. 4
pathbed=args.vcf+"/chr"+str(chr)+".vcf.gz"

#IBD status matrix
pathibd=args.ibdstatus+"/trim"+str(trim)+"_ibd2_"+str(chr)+".csv.gz"

#output file
fileibd2=args.out+"/trim"+str(trim)+"_"+"ibdal2_"+str(chr)+".csv.gz"

#Path to a bim file with the SNPs we want to include
pathbimfile=args.bim+"/chr"+str(chr)+".bim"

print ('IBD2 status file',pathibd)
print ('IBD2 output file', fileibd2)
print ('VCF file', pathbed)
print ('Chromsome ',str(chr))
print ('Bimfile ',pathbimfile)

#Read in information about ibd regions to have the siblings in the right order
#Just need the header
E = pd.read_csv(pathibd,compression='gzip',sep=',',nrows=0)
sib=E.columns
sib=sib[6:sib.shape[0]]

G=pd.read_csv(pathbed,delimiter='\t',skiprows=vcfheader,compression='gzip')
#Make sure that the columns SNP name and CHR are named that in the dataframe
G.rename(columns={'ID': 'SNP', '#CHROM': 'CHR'}, inplace=True)
#Quick fix in case the SNP column is on the format CHR,SNP
ccc=","+str(chr)
G['SNP']=G['SNP'].str.replace(ccc,"")

Eibd2=G[['SNP','REF','ALT']]

bim = pd.read_csv(pathbimfile,sep="\t",header=None)
bim.columns=["chr","SNP","b","pos","A1","A2"]
bim=bim.loc[bim['chr'].astype(str)==str(chr)]

G=G.merge(bim[['SNP']],left_on='SNP',right_on='SNP',how='inner')
G=G.sort_values(by=['POS'])

for i in range(0,len(sib)):
        s=sib[i]
        s1=s[:7]
        s2=s[8:]
        S=G[['SNP','POS',s1,s2]]
        s1ibd2=s1+"_ibd2s_"+s2
        s2ibd2=s2+"_ibd2s_"+s1
        S[s1ibd2]=(S[s1].str[0].astype(int)+S[s1].str[2].astype(int)+S[s2].str[0].astype(int)+S[s2].str[2].astype(int))*0.5
        S[s2ibd2]=S[s1ibd2]
        Eibd2=Eibd2.merge(S[['SNP',s1ibd2,s2ibd2]],left_on='SNP',right_on='SNP',how='outer')

###Write out a compressed file
Eibd2.to_csv(fileibd2,compression='gzip',index = None, header=True)
