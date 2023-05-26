# Investigating the relationship between the PGS and various phenotypes
These scripts show how we performed linear and logistic regression to attain regression results shown in *Table 1*, *Figure 3* and *Supplementary Table 1*.

## Input files
To run these scripts, the following input files are needed:
1) **GWAS summary statistics**: <br />

The T-statistics computed with the script in [/SumStatAnalysis](https://github.com/stefaniabe/PrimaryParticipationGWAS/tree/main/SumStatAnalysis).
There should be one input file per weight (BSPC, WSPC, TNTC) and the files should be tab seperated with the following header:
```
MarkerName A1 Beta
```
A1 refers to the allele that corresponds to a positive T/Z-statistic and Beta is the absolute value of the T/Z-statistic. 
In the current study, to compute the PGSs, we only included 500,632 common biallelic sequence variants (MAF>1%) that had fulfilled various quality control criteria (see Method section in paper).

2) **Genotype data**: <br />

Bed file with the genotypes of the individuals to compute PGS for. Should include all the SNPs that one wants to include in the computation of the PGS.

3) **Quantitative phenotype data**: <br />

Tab seperated matrix with quantitative phenotype values. Unique identifier column should be called "IID". The columns of the phenotypes should bear the name of each corresponding phenotype.
```
IID phenotype1 phenotype2 phenotype3
```
Values that are missing because they are above/below the reportable range should be coded as -8888888/-9999999 in the matrix. 
In the current study, the trait educational attainment had been mapped to years of schooling using the ISCED classification.

4) **Categorical phenotype data**: <br />

For categorical phenotypes, seperate files for each phenotype are needed. 
These should be tab seperated files called "phenotype.txt" and without a header. 
Unique identifier should be in the first column and the binary phenotype value (0 or 1) in the second column.

5) **Individuals to include:** <br />

Tab seperated file with unique identifiers of individuals to include in the analysis in a column with the header "IID".

6) **Covariate data:** <br />

Tab-seperated matrix with information about IID, sex, Age, YOB, PC1-PC40 and chip-type.
```
IID sex Age YOB chip PC1 PC2 ... PC40
```


7) **Phenotype names:** <br />

Various files with a list of the phenotype names to include in various types of analysis (see below).

8) **LD score regression intercepts:** <br />

Files with LD score regression intercepts attained from performing LD score regression using GWAS summary statistics based on the same set of 
individuals that we include in the regression analysis (see [/GenCorrH2](https://github.com/stefaniabe/PrimaryParticipationGWAS/tree/main/GenCorrH2)). 
The files should be called L2_all.txt (LD score regression intercepts based on GWAS of all individuals), 
L2_women.txt (LD score regression intercepts based on GWAS of women) and L2_men.txt (LD score regression intercepts based on GWAS of men). 
These should be tab seperated files and have the header "pheno" (for phenotype) and "L2" (for LD score regression intercept).
```
pheno L2
```
If one does not want to adjust the P-values with LD score regression intercepts, then all the L2 values should be set to 1 in these files.

## Computing the polygenic scores
Polygenic scores were computed with with plink [1.9](https://www.cog-genomics.org/plink/1.9/) . The inputs are as follows and in the following order (order matters):
1) bed: Path to the bed file of the individuals we want to compute PGS for.
2) gwassumstat: full path to the gwas summary statistic file.
3) pgspath: path to the folder where to store the polygenic scores.
4) pgs: Name of the pgs.

```
bash PGS.sh $bed $gwassumstat $pgspath $pgs
```
Note that the current study used the UKBB phased haplotype data and hence there were no missing genotypes in the bed file.

## Transforming quantitative phenotypes
Quantitative phenotypes were adjusted for age, YOB, sex and 40 PCs. 
The inputs for the script are as follows and in the following order (order matters):
1) phenomatrixpath: Full path to the phenotype matrix (tab seperated matrix with a header). Unique identifier column should be called "IID".
2) transpath: Path to the folder where the transformed phenotypes should be stored. Preferrably this is the same folder that the categorical phenotypes are stored.
3) indpath: Full path to a file with the unique identfiers of the individuals that we want to include in the regression analysis (header with IID). Current study included 272,409 individuals that are not closely related to each other.
4) pcpath: Full path to the covariate matrix, (header: IID Age, Sex, YOB, chip, PC1-PC40).
5) phenopath1: Full path to a tab seperated file with column called "pheno" which contains the names of the phenotypes we want to transform.
6) phenopathage: Full path to a tab seperated file with column called "pheno" which contains the names of the phenotypes we want to RINT and adjust for age up to the order to 3, sex and YOB (should have a file with just the header "pheno" if no phenotypes are to be transformed in this manner).
7) phenopathnorank: Full path to a tab seperated file with column called "pheno" which contains the names of the phenotypes we do not want to RINT but standardise, adjust for YOB up to the order of 3, sex and age (should have a file with just the header "pheno" if no phenotypes are to be transformed in this manner).
8) phenopathyob: Full path to a tab seperated file with column called "pheno" which contains the names of the phenotypes we want to RINT and adjust for YOB up to the order to 3, sex and age (should have a file with just the header "pheno" if no phenotypes are to be transformed in this manner).

```
Rscript TransformingPhenoQuant.R $phenomatrixpath $transpath $indpath $pcpath $phenopath1 $phenopathage $phenopathnorank $phenopathyob
```

## Regression analysis
Investigating the relationship between the PGS and various phenotypes (*Table 1*, *Supplementary Table 1* and *Figure 3*).
For the logistic/linear regression the inputs are as follows and in the following order:
1) transpath: Path to a folder which contains the phenotype lists, stored as individual files called "phenotype.txt". The files should be tab seperated, with no header and have the unique identifier in the first column and the phenotype value in the second column. For the quantitave phenotypes, the phenotypic values have been adjusted for Age, YOB, sex and PCs (see above).
2) outpath: Path to the folder to store the output regression tables.
3) pgsscorepath: Full path to the polygenic score file.
4) pgs: Name of the polygenic score, for naming the output files.
5) indpath: Full path to the file with the individuals that we want to include in the regression analysis (header with IID).
6) ccphenopath/phenopath: List of phenotypes we want to perform logistic/linear regression for (no header). Should correspond to the names in folder $transpath, excluding the ".txt" extension.
7) pcpath: Full path to the covariate matrix, (IID, Age, Sex, YOB, chip, PC1-PC40).
8) L2path: Folder with LDSC intercept files attained from performing LD score regression using GWAS summary statistics based on the same set of individuals ($indpath). The files should be called L2_all.txt (based on all individuals), L2_women.txt (sex-specific GWAS based on only the women) and L2_men.txt (sex-specific GWAS based on only the men) and have the header "pheno" and "L2".
9) edupath: Full path to a phenotype file with adjusted and standardised educational attainment values (see "Transforming quantitative phenotypes above").

```
Rscript Logisticregression.R $transpath $outpath $pgsscorepath $pgs $indpath $ccphenopath $pcpath $L2path $edupath
Rscript Linearregression.R $transpath $outpath $pgsscorepath $pgs $indpath $phenopath $pcpath $L2path $edupath
```

## Making figure
Making a figure that shows regression results from investigating the relationship between the PGS and various phenotypes in the subset of men and women seperately (*Figure 3*). 
The inputs are as follows and in the following order:
1) pgs: Name of the pgs for the output
2) WL: Full path to the file with the women specific linear regression output from above
3) WC: Full path to the file with the women specific logistic regression output from above
4) ML: Full path to the file with the men specific linear regression output from above
5) MC: Full path to the file with the men specific logistic regression output from above
6) figpath: Path to the folder to store the figure

```
Rscript SexDiffFig.R $pgs $WL $WC $ML $MC $figpath
```
