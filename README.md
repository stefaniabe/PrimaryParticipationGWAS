# Reproducing the main analyses in *Studying the genetics of participation using footprints left on the ascertained genotypes*
Scripts to reproduce the main analyses in Benonisdottir and Kong (2022) (https://doi.org/10.1101/2022.02.11.480067). <br />
* **/GenCorrH2**: Performing LD score regression with the program [LDSC version 1.0.1](https://github.com/bulik/ldsc) 
to attain LD score regression intercepts, genetic correlation estimates and inputs for heritability estimates on a liability scale - *Table 2*, *Supplementary Table 2* and *Supplementary Table 3*.
* **/PartGWAS**: Performing primary participation GWAS. <br />
* **/PGSPheno**: Performing linear and logistic regression to examine the relationship between the primary participation polygenic score 
and various phenotypes. - *Table1*, *Figure 3* and *Supplementary Table 1*.
* **/Simulations**: Simulations to investigate the relative allele frequencies in different groups of individuals (participants, non-participants, individuals with a participating/non-participating sibling) 
and genetic segments (IBD shared and not-shared) - *Figure 4* and *Extended Data Figure 8*. 
* **/SumStatAnalysis**: Computing T-statistics from the summary statistics files generated with scripts in /PartGWAS and investigating 
their relationship with allele frequency - *Extended Data Figure 5*, *Extended Data Figure 6* and *Extended Data Figure 7*. <br />
