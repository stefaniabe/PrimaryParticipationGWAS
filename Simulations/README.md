# Simulations
This script shows how we performed simulations in R (version 3.4.3) assuming a liability-threshold model to investigate, for a given SNP, allele frequencies in different groups of individuals 
(participants, non-participants, individuals with a participating/non-participating sibling) and genetic segments - *Figure 4* and *Extended Data Figure 8*.
The input parameters are as follows and in the following order:
1) N: Number of individuals in the population
2) a: Target participation rate
3) outpath: Path to folder to store outputs (table and figure)
4) name: Short name for the output matrix
5) p: Allele frequency of the SNP in question
6) sim: Number of simulations

```
Rscript LiabilityModel.R $N $a $outpath $name $p $sim
```

Note in the paper the simulations were repeated 500 times before the figures were made.
If computational resources are limited one might want to simulate first many output matrices parallel (choosing a different $name for each run), join them in the end and then make the figures.
