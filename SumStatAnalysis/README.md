# Summary Statistic Analysis
This script shows how we used the output summary statistics files attained 
with scripts in [/partGWAS](https://github.com/stefaniabe/PrimaryParticipationGWAS/tree/main/PartGWAS) to  <br />
a) Compute BSPC, WSPC and TNTC T-statistics.  <br />
b) Adjust the WSPC and TNTC T-statistics with a two-step step allele frequency adjustment.  <br />
c) Investigate relationship with allele frequency - *Extended Data Figures 5, 6 and 7*.  <br />

The input parameters are as follows and in the following order:
1) "trim250_AF_sumstat.csv.gz": Full path to the summary statistic file showing sibling IBD specific summary statistics. Current study used the trim250 analysis.
2) "TNT_AF_sumstat.csv.gz": Full path to the summary statistic file showing parent-offspring transmitted and non-transmitted summary statistics.
3) "TNT_DH_sumstat.csv.gz": Full path to the summary statistic file that shows the count of errors of inferring the shared allele for double-heterosygous parent-offspring pairs.
4) "hiqhqualitySNPs.txt": Full path to a tab seperated file that has a column called "SNP" which contains the SNP IDs of common high-quality SNPs. Current study used the 500,632 SNPs used to compute the primary participation PGS.
5) "outpath": Path to the folder where the outputs should be stored

```
Rscript Tstatistics.R "trim250_AF_sumstat.csv.gz" "TNT_AF_sumstat.csv.gz" "TNT_DH_sumstat.csv.gz" "hiqhqualitySNPs.txt" "outpath"
```
