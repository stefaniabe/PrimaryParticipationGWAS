# LD score regression: intercepts, genetic correlations and heritability
These scripts show how we used the program [LDSC](https://github.com/bulik/ldsc) version 1.0.1 to attain:  
A) LD score regression intercepts used to adjust P-values shown in *Table 1*. <br />
B) Genetic correlations shown in *Table 2* and *Supplementary Table 2*. <br />
C) Heritablity estimates for categorical traits that were then transformed to liabilty score heritability estimates (see Methods in paper) and shown in *Supplementary Table 3*. <br />

The LSR analysis included the following input files:
1) **GWAS summary statistics**: <br />

These are the BSPC and WSPC T-statistics (see [/SumStatAnalysis](https://github.com/stefaniabe/PrimaryParticipationGWAS/tree/main/SumStatAnalysis)) and GWAS summary statistics from performing 'traditional' GWAS for the traits shown in Table 1. 
The 'traditional' GWAS results were attained with [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/) and the P_LINREG column in the output was used as input P-value to LDSC. <br />

The required format of the input GWAS summary statistics files is well described on the [wiki](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format) page of the LDSC program.
In the current study, the input GWAS summmary statistics had been merged with the 500,632 SNPs that were used to compute the primary participation PGS and converted to the LDSC format using munge_sumstats.py. The resulting file had the following 5 columns:
```
SNP  A1  A2  Z  N
```
with A1 being the effect allele, Z being the signed Z-score corresponding to the association and N being the sample size.
For categorical traits, N was equal to cases+controls. For the BSPC statistics the sample size for each SNP was the number of IBD2 pairs plus twice the number of IBD0 pairs. 
For the WSPC statistics N was equal to 1.5 times the number of IBD1 pairs. <br />

2) **LD scores**: <br />

For the current study, UKBB EUR LD scores were downloaded from https://pan.ukbb.broadinstitute.org on April 27th 2021 (Pan-UKB team, 2020). The required format of the input LD score files is well described on the [wiki]([https://github.com/bulik/ldsc/wiki/](https://github.com/bulik/ldsc/wiki/LD-File-Formats)) page of the LDSC program. 
The LD score input file was a zipped file ending with ".l2.ldscore.gz" and had the following columns:
```
CHR SNP BP CM MAF L2
```
In the same folder was another file with the same prefix but ending with ".l2.M_5_50" which had information about the number of variants with MAF>0.05.

## Genetic correlations
Genetic correlations were estimated with a constrained intercept LD score regression. 
That is, the intercept for the covariance was set to zero as there was no overlap of individuals between the different studies. 

The inputs are as follows and in the following order:
1) insumstat1: full path to the gwas summary statistic file of trait 1 (BSPC or WSPC).
2) insumstat2: full path to the gwas summary statistic file of trait 2 (the 'traditional' GWAS results)
3) ldpath: path to the folder and prefix of the LD score files (omitting ".l2.ldscore.gz" and ".l2.M_5_50").
4) outpath: path to where to store the output, folder and prefix

```
bash gencor.sh $insumstat1 $insumstat2 $ldpath $outpath
```
## Heritability and LD score regression intercepts
Heritablity was estimated for BSPC and the secondary participation traits with LDSC and 
then the output was converted to a liability scale as is described in the Method section of the paper.

The inputs are as follows and in the following order:
1) insumstat1: full path to the gwas summary statistic file of the trait
2) ldpath: path to the folder and prefix of the LD score files (omitting ".l2.ldscore.gz" and ".l2.M_5_50").
3) outpath: path to where to store the output, folder and prefix

For BSPC we used a constrained intercept regression (the intercept was constrained to 1):
```
bash h2_constrain.sh $insumstat1 $ldpath $outpath
```
For the other traits we did not constrain the intercept:
```
bash h2.sh $insumstat1 $ldpath $outpath
```
For the binary phenotypes, outputs from these analysis were then converted to a liability scale heritability estimate as is described in the Method section of the paper: "Estimating heritability for liability scores".
The unconstrained analysis was performed for all the traits in Table 1 in order to attain estimateds of the LD score regression intercepts.
