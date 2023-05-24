import argparse
import pandas as pd
import numpy as np
from pandas import DataFrame

argParser = argparse.ArgumentParser()
argParser.add_argument("--out", help="Path to the folder where the outputfile will be stored.")
argParser.add_argument("--ibdmatrix", help="Path to the folder where the input IBD genotype matrix (output from a previous step) is stored.")
argParser.add_argument("--trios", help="#Full path to file with trios to include (columns ID1 and ID2, each offspring will appear twice in column ID1 as these are trios).")
argParser.add_argument("--chr", help="Chromosome")

args = argParser.parse_args()
print("args=%s" % args)

chr=args.chr
#path to input
path1=args.ibdmatrix
#Path to folder to store output
path3=args.out
#Full path to file with trios to include (columns ID1 and ID2, each offspring will appear twice in column ID1 as these are trios).
pairs=args.trios

chrc="chr"+str(chr)
ibdfn=path1+"/TNT_al1n_"+str(chr)+".csv.gz"
ibdfs=path1+"/TNT_al1s_"+str(chr)+".csv.gz"
fibd=path3+"/phase_error_trios_"+str(chr)+".csv.gz"

##Need to make shared and not shared parent-offspring
pair=pd.read_csv(pairs,delimiter="\t")
pair1s=pair['ID1'].astype(str)+"_ibd1s_"+pair['ID2'].astype(str)
pair1n=pair['ID1'].astype(str)+"_ibd1n_"+pair['ID2'].astype(str)

Gibn=pd.read_csv(ibdfn,compression='gzip')
Gibs=pd.read_csv(ibdfs,compression='gzip')

#offspring
oibn=Gibn[list(pair1n)].copy()
oibs=Gibs[list(pair1s)].copy()

###Order all the matrices in a alphabetical number (the columns)
oibn=oibn.sort_index(axis=1)
oibs=oibs.sort_index(axis=1)

#Make sure that the columns are in same order for the parents
onames=oibs.columns
pps=[None]*len(onames)
ppn=[None]*len(onames)

for i in range(0,len(onames)):
        s1=onames[i]
        pp1=s1.split("_ibd1s_")[1]
        pp2=s1.split("_ibd1s_")[0]
        pps[i]=pp1+"_ibd1s_"+pp2
        ppn[i]=pp1+"_ibd1n_"+pp2

#parent
pibn=Gibn[list(ppn)].copy()
pibs=Gibs[list(pps)].copy()

#Discordant genotypes were defined as 0.5, do not want to include them in this analysis.
pibn[pibn==0.5]=777
pibs[pibs==0.5]=777
oibn[oibn==0.5]=777
oibs[oibs==0.5]=777

#Want to have three matrices, offspring, parent 1 and parent 2
#to know where we have triple heterozygous, double heterozygous and no heterozygous
o1s=oibs[oibs.columns[::2]]
o2s=oibs[oibs.columns[1::2]]
o1n=oibn[oibn.columns[::2]]
o2n=oibn[oibn.columns[1::2]]

p1s=pibs[pibs.columns[::2]]
p2s=pibs[pibs.columns[1::2]]
p1n=pibn[pibn.columns[::2]]
p2n=pibn[pibn.columns[1::2]]

#parent genotype
P1=pd.DataFrame.as_matrix(p1s)+pd.DataFrame.as_matrix(p1n)
P2=pd.DataFrame.as_matrix(p2s)+pd.DataFrame.as_matrix(p2n)

#offspring genotype
G=pd.DataFrame.as_matrix(o1s)+pd.DataFrame.as_matrix(o1n)

#When are parent or offspring heterozygous
P1[P1!=1]=0
P2[P2!=1]=0
G[G!=1]=0

#When are parent and offspring both heterozygous
PG1=P1+G
PG1[PG1!=2]=0
PG1[PG1==2]=1

PG2=P2+G
PG2[PG2!=2]=0
PG2[PG2==2]=1

#When is offspring and one parent heterozygous but not the other parent
PG=PG1+PG2
##Set triple heterozygous cases at zero
PG[PG==2]=0
#count trios for which offspring is heterozygous, one parent is heterozygous and the other parent is homozygous
ntrios=PG.sum(axis=1)

#Number of dh PO pairs for which the other parent is homozygous 00 or 11 (want to differentiate)
#parent 1 genotype
P1=pd.DataFrame.as_matrix(p1s)+pd.DataFrame.as_matrix(p1n)
P1=P1+700
#Multiply with a matrix which shows when the offspring and other parent are both heterozygous
p1=np.multiply(PG2,P1)
p1=pd.DataFrame(p1)
p1[p1==0]=np.nan
p1=p1-700
p1[p1>700]=np.nan

pa1=p1.copy()
pa1[pa1==1]=np.nan
pa1[pa1==0]=np.nan
pa1[pa1==2]=1
nh11=pa1.sum(axis=1)

pb1=p1.copy()
pb1[pb1==1]=np.nan
pb1[pb1==2]=np.nan
pb1[pb1==0]=1
nh01=pb1.sum(axis=1)

#parent 2 genotype
P2=pd.DataFrame.as_matrix(p2s)+pd.DataFrame.as_matrix(p2n)
P2=P2+700
p2=np.multiply(PG1,P2)
p2=pd.DataFrame(p2)
p2[p2==0]=np.nan
p2=p2-700
p2[p2>700]=np.nan

pa2=p2.copy()
pa2[pa2==1]=np.nan
pa2[pa2==0]=np.nan
pa2[pa2==2]=1
nh12=pa2.sum(axis=1)

pb2=p2.copy()
pb2[pb2==1]=np.nan
pb2[pb2==2]=np.nan
pb2[pb2==0]=1
nh02=pb2.sum(axis=1)

#combined
#other parent is homozygous 1
nh1=nh11+nh12
#other parent is homozyogus 0
nh0=nh01+nh02

#When one parent is homozygous, shared alle for the double heterozygous pairs can be determined withoug phasing information
os=pd.DataFrame.as_matrix(o1s)-pd.DataFrame.as_matrix(o2n)
os=os+700
Os=np.multiply(PG,os)
Os=pd.DataFrame(Os)
Os[Os==0]=np.nan
Os=Os-700

#Number of errors when the other parent is homozygous 1
Osn1=Os.copy()
Osn1[Osn1!=1]=0
ne1=Osn1.sum(axis=1)

#Number of errors when the other parent is homozygous 0
Osn0=Os.copy()
Osn0[Osn0!=-1]=0
Osn0[Osn0==-1]=1
ne0=Osn0.sum(axis=1)

#To make the sumstat matrix
G=Gibn[['SNP','REF','ALT']].copy()
G['n_DH_other1']=nh1
G['n_DH_other0']=nh0
G['error_1']=ne1
G['error_0']=ne0

G.to_csv(fibd,compression ='gzip',index = None, header=True)
