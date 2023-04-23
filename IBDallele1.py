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

chr=args.chr
vcfheader=int(args.headerline)-1

#This script makes an output matrix for the sibling genotypes as if all of SNPs are shared IBD1
#In the next step, this output matrix is multiplied with the IBD status matrix to produce a matrix 
#where we can see genotypes that are truly shared IBD1 
#other values will be NA.
#This is run for trim=0, as the trimming step is done in the IBDstatus step

#Run this with trim=0. 
#Only need to run this script once no matter the trim, as trimming information is in the IBD status matrix.
#This saves future computing time if one wants to run this for more trimming scenarios
trim=0

#Path to ibd1 status matrix
pathibd=args.ibdstatus+"/trim"+str(trim)+"_ibd1_"+str(chr)+".csv.gz"

#File with chromosome specific vcf files
#The vcf file has the header in row nr. 4
pathbed=args.vcf+"/chr"+str(chr)+".vcf.gz"

#Path to a bim file with the SNPs we want to include
pathbimfile=args.bim+"/chr"+str(chr)+".bim"

#The output files
fileibd1s=args.out+"/trim"+str(trim)+"_ibdal1s_"+str(chr)+".csv.gz"
fileibd1n=args.out+"/trim"+str(trim)+"_ibdal1n_"+str(chr)+".csv.gz"

print ('IBD1 status file',pathibd)
print ('IBD1 shared output file', fileibd1s)
print ('IBD1 not-shared output file', fileibd1n)
print ('VCF file', pathbed)
print ('Chromsome ',str(chr))
print ('Bimfile ',pathbimfile)

####Read in information about ibd regions
E = pd.read_csv(pathibd,compression='gzip',sep=',')
sib=E.columns
#First 6 columns have information about the SNP in question, the other columns have sibling-pair IDs
sib=sib[6:sib.shape[0]]

G=pd.read_csv(pathbed,delimiter='\t',skiprows=vcfheader,compression='gzip')
#Make sure that the columns SNP name and CHR are named that in the dataframe
G.rename(columns={'ID': 'SNP', '#CHROM': 'CHR'}, inplace=True)
#Quick fix in case the SNP column is on the format CHR,SNP
ccc=","+str(chr)
G['SNP']=G['SNP'].str.replace(ccc,"")

bim = pd.read_csv(pathbimfile,sep="\t",header=None)
bim.columns=["chr","SNP","b","pos","A1","A2"]
bim=bim.loc[bim['chr'].astype(str)==str(chr)]

G=G.merge(bim[['SNP']],left_on='SNP',right_on='SNP',how='inner')
G=G.sort_values(by=['POS'])

Eibd1s=G[['SNP','REF','ALT']]
Eibd1n=G[['SNP','REF','ALT']]

for i in range(0,len(sib)):
        s=sib[i]
        s1=s[:7]
        s2=s[8:]
        S=G[['SNP','POS',s1,s2]].merge(E[['SNP',s]],left_on='SNP',right_on='SNP')
        s1ibd1s=s1+"_ibd1s_"+s2
        s2ibd1s=s2+"_ibd1s_"+s1
        s1ibd1n=s1+"_ibd1n_"+s2
        s2ibd1n=s2+"_ibd1n_"+s1
        ##Which alleles are shared and which are not
        ST=S.loc[S[s]==1]
        if not ST.empty:
                ST1=ST.loc[(ST[s1].str[0].astype(int)+ST[s1].str[2].astype(int)+ST[s2].str[0].astype(int)+ST[s2].str[2].astype(int))==3]
                ST1['a1']=np.where(ST1[s1].str[:3]=="1|0",555,np.where(ST1[s1].str[:3]=="0|1",666,np.nan))
                ST1['a2']=np.where(ST1[s2].str[:3]=="1|0",555,np.where(ST1[s2].str[:3]=="0|1",666,np.nan))
                ST2=ST.loc[(ST[s1].str[0].astype(int)+ST[s1].str[2].astype(int)+ST[s2].str[0].astype(int)+ST[s2].str[2].astype(int))==1]
                ST2['a1']=np.where(ST2[s1].str[:3]=="1|0",666,np.where(ST2[s1].str[:3]=="0|1",555,np.nan))
                ST2['a2']=np.where(ST2[s2].str[:3]=="1|0",666,np.where(ST2[s2].str[:3]=="0|1",555,np.nan))
                ST = ST1.append(ST2)
                ST['POS']=ST['POS'].astype(float)
                S=S.merge(ST[['SNP','a1','a2']],left_on='SNP',right_on='SNP',how="outer")
                S['POS']=S['POS'].astype(float)
                S=S.sort_values(by=['POS'])
                ##Need to find the first and last non-na values so that I can have them at the top and bottom of the column
                ##Otherwise interpolate will not fill in the values like it should
                a1nf=S['a1'].first_valid_index()
                a1nl=S['a1'].last_valid_index()
                a2nf=S['a2'].first_valid_index()
                a2nl=S['a2'].last_valid_index()
                a1nfs=str(S['a1'].first_valid_index())
                a2nfs=str(S['a2'].first_valid_index())
                ##Size of matrix
                m=S.shape[0]-1
                if a1nfs!="None":
                        S['a1'].loc[0]=S['a1'].loc[a1nf]
                        S['a1'].loc[m]=S['a1'].loc[a1nl]
                else:
                        S['a1'].loc[0]=777
                        S['a1'].loc[m]=777
                if a2nfs!="None":
                        S['a2'].loc[0]=S['a2'].loc[a2nf]
                        S['a2'].loc[m]=S['a2'].loc[a2nl]
                else:
                        S['a2'].loc[0]=777
                        S['a2'].loc[m]=777
                ##Fill in with interpolate
                S['a2']=S['a2'].interpolate(method='nearest',limit=None)
                S['a1']=S['a1'].interpolate(method='nearest',limit=None)
                S['b']=S[s1].str[0].astype(int)+S[s1].str[2].astype(int)+S[s2].str[0].astype(int)+S[s2].str[2].astype(int)
                S['c']=abs(S[s1].str[0].astype(int)+S[s1].str[2].astype(int)-S[s2].str[0].astype(int)-S[s2].str[2].astype(int))
                S[s1ibd1s]=np.where((S['b'])>2,1,
                np.where(S['b']<2,0,
                np.where(S['c']==2,0.5,
                np.where(S['a1']==555,S[s1].str[0].astype(int),
                np.where(S['a1']==666,S[s1].str[2].astype(int),
                np.where(S['a2']==555,S[s2].str[0].astype(int),
                np.where(S['a2']==666,S[s2].str[2].astype(int),888)))))))
                S[s2ibd1s]=S[s1ibd1s]
                S['d1']=S[s1].str[0].astype(int)+S[s1].str[2].astype(int)-S[s1ibd1s]
                S['d2']=S[s2].str[0].astype(int)+S[s2].str[2].astype(int)-S[s2ibd1s]
                S[s1ibd1n]=np.where(S['c']==2,0.5,S['d1'])
                S[s2ibd1n]=np.where(S['c']==2,0.5,S['d2'])
        else:
                #If a pair has no IBD=1 region at this chromosome, set all the values at 5, will be "NA" when multiplied with the IBD1status matrix in the end
                S[s1ibd1n]=5
                S[s1ibd1s]=5
                S[s2ibd1n]=5
                S[s2ibd1s]=5
        Eibd1s=Eibd1s.merge(S[['SNP',s1ibd1s,s2ibd1s]],left_on='SNP',right_on='SNP',how='outer')
        Eibd1n=Eibd1n.merge(S[['SNP',s1ibd1n,s2ibd1n]],left_on='SNP',right_on='SNP',how='outer')

###Write out compressed files
Eibd1s.to_csv(fileibd1s,compression='gzip',index = None, header=True)
Eibd1n.to_csv(fileibd1n,compression='gzip',index = None, header=True)                                                                       
