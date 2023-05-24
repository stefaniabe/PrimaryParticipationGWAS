# PrimaryParticipationGWAS
Python and R scripts to perform primary participation GWAS as in Benonisdottir and Kong (2022) (https://doi.org/10.1101/2022.02.11.480067).
This is done by comparing IBD shared and not-shared alleles of participating first-degree relatives and is performed separately for sibling pairs and parent-offspring pairs. 

# Input files and requirements

The python scripts require python 2.7, argparse, numpy and pandas.

To perform the analyses the following datafiles are required:

1) Phased genotypes of the first-degree relatives on a vcf format. These should be chromosome specific, tab separated files named chr*.vcf.gz. The header should be as follows:
```
#CHROM  POS   ID REF ALT QUAL FILTER INFO FORMAT        ID1           ID2      ...    IDN
```
The header does not necessarily need to be in the first line. However, the line number of the header should be specified as an input variable when the scripts requiring the vcf files are run. 
E.g. when converting the phased bgen files available in the UKBB data release 2 to vcf format with QCTOOL 2.0.1 (see https://www.well.ox.ac.uk/~gav/qctool_v2/), one gets a vcf file with the header in the fourth line and hence the headerline input should be 4. Columns ID1-IDN are the unique identifiers of the N individuals that constitue the first degree relatives. ALT stands for the allele coded as 1. The SNP IDs (the column called ID) need to be unique per chromosome.
For the individual level columns (ID1-IDN), only the first three values in each field are used, and they should be either 0|0, 0|1, 1|0 or 1|1. E.g. a line in the input file could be as follows:
```
 1   7527219 SNP1 G   A   .     .     .    GT:GP    0|1,1,0,0,1  0|1,1,0,0,1   ... 1|1,0,1,0,1
```
but the script would interpret it as:
```
 1   7527219 SNP1 G   A   .     .     .    GT:GP       0|1            0|1     ...     1|1 
```

2) Bimfiles with the SNPs that are in the vcf files. These should be chromosome specific tab seperated files named chr*.bim with no header but the following columns:
```
CHR SNP CM BP A1 A2
```
Given that the purpose of this file is to ensure that we are always working with the same SNPs throughout the steps and to order the SNPs by BP position, CM column can be populated with zeros. The SNP IDs (the SNP column) need to be unique per chromosome.

## Sibling pair specific files
3) The IBD segments of the sibling pairs as inferred by the program snipar (https://github.com/AlexTISYoung/snipar, see Young et. al (2022) https://doi.org/10.1038/s41588-022-01085-0).
These should be chromosome specific, tab separated files which are named chr_*.ibd.segments.gz and have the following header:
```
ID1 ID2 IBDType Chr start_coordinate stop_coordinate startSNP stopSNP length
```

4) A file with the unique identifiers of the sibling pairs that we want to compute IBD specific allele frequencies for. Note, this can be a subset of the siblings in the vcf file mentioned above and enables you to compute IBD specific allele frequencies for different subgroups of sibling pairs.
This file should be a tab seperated and have the following header:
```
fid ID1 ID2
```
fid stands for family id and ID1 and ID2 are the unique personal identifiers for the two siblings (one row per pair).
It should also be noted that the variances computed in step 4 in the sibling protocol assume that the shared versus not-shared frequency differences of each pair are independent of each other and that for every SNP, the allele frequencies computed for individual sibling pairs, whether IBD0 or IBD2 pairs are independent of each other. In the paper, we chose the two first participating siblings in sib-ships that were larger than two, and hence the same individual was only observed once in each pair.

5) Minor allele frequency (MAF) files containing information about MAF computed in the group of the sibling pairs in file 4 for the SNPs that are in the vcf and bim files. These should be tab separated, chromsome specific files named chr*.frq with the following header:
```
CHR SNP A1 A2 MAF NCHROBs
```
Note, the A1 allele is the minor allele.

## Parent-offspring specific files

6) A file with the unique identifiers of the parent-offspring pairs that we want to compute transmitted and non-transmitted allele frequencies for.
This file should be tab separated, and have the following header:
```
fid ID1 ID2
```
fid stands for family id. ID1 is the identifier for the offspring and ID2 is the identifier for the parent. It should be noted that the variances computed in step 3 in the parent-offspring protocol assume that the transmitted versus non-transmitted frequency differences of each pair are independent of each other. In the paper, we chose to include parent-offspring pairs for which the offspring did not have a participating sibling, and hence the same individual was only observed once among the pairs.

7) MAF files containing information about MAF computed in the group of the parents (ID2 in file 6) for the SNPs that are in the vcf files and bim files. These should be tab separated, chromsome specific files named chr*.frq with the following header:
```
CHR SNP A1 A2 MAF NCHROBs
```
Note, the A1 allele is the minor allele.

# Sibling pairs: Stepwise protocol
The stepwise protocol listed below serves to compute IBD specific allele frequencies for sibling pairs. It assumes that sibling-pairs have been identified and that their IBD shared segments have been inferred with the program snipar. To run all the scripts, datafiles 1)-5), in the section 'Input files and requirements' above, are required.

## Step 1: IBD status

For sibling pairs the IBD status scripts ( IBDstatus0.py, IBDstatus1.py and IBDstatus2.py) should be run to compute matrices with information about which SNPs each sibling pair shares IBD0, IBD1 and IBD2. Each of those scripts produces an identity matrix with one row per SNP and columns corresponding to sibling pairs. For each IBD matrix, a value of 1 indicates that a given sibling pair shares that SNP of the correspondent IBD status, with a value of zero indicating the opposite.

To account for uncertainty in the IBD inferring algorithm close to recombination events, SNPs can be trimmed from the start and end of each IBD region for each sibling pair (input --trim). For subsequent analysis it is necessary that trim is set to 0 (--trim 0) in step 1 in the first round, but it is recommended that you try out other trimming scenarios as well, e.g. --trim 100, --trim 250 and --trim 500, and compare those results. Step 3, 4 and 5 also require --trim to be specified.

The following 5 arguments are required:
1) --chr: Chromosome
2) --ibd: Path to the folder where with the IBD-inferred segments as inferred by snipar are stored (chromosome specific compressed files with the name chr_*.ibd.segments.gz)
3) --bim: Path to the folder where the bim-files with the SNPs in the vcf files are stored (the chromosome specific files should be named chr*.bim)
4) --trim: Number of SNPs to trim from the beginning and end of the IBD regions.
5) --out: Path to the folder where the output will be stored.

Example:
```
python IBDstatus0.py --chr 22 --ibd "ibdsegpath" --bim "bimpath" --trim 0 --out "ibdstatuspath"
python IBDstatus1.py --chr 22 --ibd "ibdsegpath" --bim "bimpath" --trim 0 --out "ibdstatuspath"
python IBDstatus2.py --chr 22 --ibd "ibdsegpath" --bim "bimpath" --trim 0 --out "ibdstatuspath"
```

## Step 2: IBD shared allele
 
For sibling pairs, the IBDallele scripts, IBD0allele.py, IBD1allele.py and IBD2allele.py, make output matrices for the sibling-pair genotypes under the assumption that all of SNPs are shared IBD0, IBD1 and IBD2 respectively (which they are not). In step 3, these IBDallele matrices are multiplied with the IBDstatus matrices generated in step 1. Please note that the output files from this step are only intermediate files that are required as an input in step 3. They should not be used otherwise.

Of all the steps, step 2 requires the most time and resources. To save computing time, it is here assumed that trim=0, as the trimming is done in step 1. Thus, step 2 only needs to be performed once independent of how many trimming scenarios one wants to explore.

The following 5 arguments are required:
1) --chr: Chromosome
2) --ibdstatus: Path to the folder where the ibd status matrix (output from the previous step) was stored.
3) --bim: Path to the folder where the bim-files with the SNPs in the vcf files are stored (the chromosome specific files should be named chr*.bim)
4) --vcf: Path to the folder with the chromosome specific zipped VCF with the phased genotypes of individuals with a sibling in the dataset (compressed files called chr*.vcf.gz)
5) --headerline: Line number of the header (column names) in the vcf files.
6) --out: Path to the folder where the output is written to.

Example:
```
python IBDallele0.py --chr 22 --ibdstatus "ibdstatuspath" --bim "bimpath" --vcf "vcfpath" --headerline 4 --out "ibdallelepath"
python IBDallele1.py --chr 22 --ibdstatus "ibdstatuspath" --bim "bimpath" --vcf "vcfpath" --headerline 4 --out "ibdallelepath"
python IBDallele2.py --chr 22 --ibdstatus "ibdstatuspath" --bim "bimpath" --vcf "vcfpath" --headerline 4 --out "ibdallelepath"
```
 
## Step 3: IBD matrix multiplication
Here, the output matrices from IBDstatus (step 1) and IBDallele (step 2) are multiplied, to produce matrices that, for each sibling pair, show genotype information only for SNPs shared IBD0, IBD1 and IBD2 respectively. Other values in the matrices will be NA. E.g. if SNP in row i for sibling pair j is shared IBD2, then the IBD2 matrix will display the shared genotype at row i and column j, else the value will be NA.

The following 6 arguments are required:
1) --chr: Chromosome
2) --ibdstatus: Path to the folder where the input IBD status matrix (output from step 1) is stored
3) --ibdallele: Path to the folder where the input IBD allele matrix (output from step 2) is stored
4) --trim: Number of SNPs to trim from the beginning and end of the IBD regions. (This is to make sure that the ibdallele matrix is multiplied with the right ibdstatus matrix).
5) --maf: Path to the folder with the files that show MAF for the SNPs in the vcf and bim files. The MAF is only used for imputation for IBD1 when it is impossible to infer the shared and not-shared allele for double heterozygous siblings. This should happen rarely, especially with a dense set of SNPs. E.g. in the UKBB phased genotype array data (approx 660,000 SNPs), the shared and not-shared allele for IBD1 could always been inferred from phasing information for double heterozygous sibling pairs.
6) --out: Path to where you want to store the output.

Example:
```
python IBDmult0.py --chr 22 --ibdstatus "ibdstatuspath" --ibdallele "ibdallelepath" --trim 0 --out "ibdmatrixpath"
python IBDmult1n.py --chr 22 --ibdstatus "ibdstatuspath" --maf "mafpath" --ibdallele "ibdallelepath" --trim 0 --out "ibdmatrixpath"
python IBDmult1s.py --chr 22 --ibdstatus "ibdstatuspath" --maf "mafpath" --ibdallele "ibdallelepath" --trim 0 --out "ibdmatrixpath"
python IBDmult2.py --chr 22 --ibdstatus "ibdstatuspath" --ibdallele "ibdallelepath" --trim 0 --out "ibdmatrixpath"
```

## Step 4: IBD specific allele frequencies and variances
Computing IBD specific allele frequencies and variances from the output file produced in step 3.

The scripts requires the following 5 arguments:
1) --chr: Chromosome
2) --trim: Number of SNPs to trim from the beginning and end of IBD regions for each sibling pair
3) --ibdmatrix: Path to the folder where the input IBD genotype matrix (output from step 3) is stored
4) --siblings: Path to a file with the sibling pairs you want to include in your allele frequency computations. The file should be tab delimieted, have three columns and have the header "fid ID1 ID2". Each line corresponds to one sibling pair. fid: Family ID, ID1: identifier for sibling 1, ID2: identifier for sibling 2.
5) --out: Path to where you want to store the SNP wise output, IBD specific allele frequencies and variances.

Example:
```
python IBDaf0.py --chr 22 --trim 0 --ibdmatrix "ibdmatrixpath" --siblings "pairstoincludepath" --out "ibdfreqpath"
python IBDaf1n.py --chr 22 --trim 0 --ibdmatrix "ibdmatrixpath" --siblings "pairstoincludepath" --out "ibdfreqpath"
python IBDaf1s.py --chr 22 --trim 0 --ibdmatrix "ibdmatrixpath" --siblings "pairstoincludepath" --out "ibdfreqpath"
python IBDaf2.py --chr 22 --trim 0 --ibdmatrix "ibdmatrixpath" --siblings "pairstoincludepath" --out "ibdfreqpath"
python IBDvar0.py --chr 22 --trim 0 --ibdmatrix "ibdmatrixpath" --siblings "pairstoincludepath" --out "ibdvarpath"
python IBDvar1.py --chr 22 --trim 0 --ibdmatrix "ibdmatrixpath" --siblings "pairstoincludepath" --out "ibdvarpath"
python IBDvar2.py --chr 22 --trim 0 --ibdmatrix "ibdmatrixpath" --siblings "pairstoincludepath" --out "ibdvarpath"
```

## Step 5: Make the summary statistic file
Combines the output files from step 4 in one summary statistic file.
The script requires R and the following 5 arguments in this order:
1) trim: Number of SNPs to trim from the beginning and end of IBD regions for each sibling pair (this is for name purposes).
2) ibdfreqpath: Path to the folder with the IBD specific allele frequencies (output of step 4)
3) ibdvarpath: Path to the folder with the IBD specific variances (output of step 4)
4) bimpath: Path to the folder with the BIM files
5) mafpath: Path to the folder with the MAF files
6) summarystatisticpath: Path to where to store the summary statistic file

Example:
```
Rscript IBDafsumstat.R 0 "ibdfreqpath" "ibdvarpath" "bimpath" "mafpath" "summarystatisticpath"
```
## Output file
This stepwise protocol results in a summary statistic file that has a name starting with 'trim0_AF_sumstat' and has the following columns:
1) SNP: SNP ID
2) CHR: Chromosome
3) POS: Position
4) AF: Allele frequency of the ALT allele as given in the MAF file
5) REF: Allele coded as 0
6) ALT: Allele coded as 1
7) ibd0: Number of sibling pairs underlying the ib0mean computations
8) ib0mean: Allele frequency among sibling pairs sharing the SNP IBD0
9) ib0meanvar: Variance of the allele frequency per sibling pair among the sibling pairs sharing the SNP IBD0
10) ibd1n: Number of sibling pairs underlying the ib1nmean computations
11) ib1nmean: Allele frequency for the non-shared alleles among the sibling pairs sharing the SNP IBD1
12) ibd1s: Number of sibling pairs underlying the ib1smean computations (should be the same as ibd1n).
13) ib1smean: Allele frequency for the shared alleles among the sibling pairs sharing the SNP IBD1
14) ib1diffvar: Variance of the differene between the shared allele and the average non-shared allele per pair.
15) ibd2: Number of sibling pairs underlying the ib2mean computations
16) ib2mean: Allele frequency among sibling pairs sharing the SNP IBD2
17) ib2meanvar: Variance of the allele frequency per sibling pair among the sibling pairs sharing the SNP IBD2

# Parent-offspring pairs: Stepwise protocol
Here it is shown how the transmitted and non-transmitted allele frequencies are computed for parent-offspring pairs.
The stepswise protocol listed below assumes that parent-pairs have been identified.

To run all the scripts, datafiles 1), 2), 6) and 7), described above in the section "Input files", are required.

## Step 1: Transmitted and non-transmitted allele
 
Infer which alleles are transmitted and non-transmitted for the PO-pairs.

The following 7 arguments are required:
1) --chr: Chromosome
2) --bim: Path to the folder where the bim-files with the SNPs in the vcf files are stored (the chromosome specific files should be named chr*.bim)
3) --vcf: Path to the folder with the chromosome specific zipped VCF with the phased genotypes of the parent-offspring pairs (compressed files called chr*.vcf.gz)
4) --maf: Path to the folder with the files that shows minor allele frequency (MAF) for the SNPs in the vcf files. The MAF is only used for imputation when it is impossible to infer the transmitted and non-transmitted allele for double heterozygous parent-offspring pairs. This should happen rarely, especially with a dense set of SNPs.
5) --po_pairs:Path to a file with the parent-offspring pairs you want to include in your allele frequency computations. The file should be tab delimieted, have three columns and have the header "fid ID1 ID2". Each line corresponds to one parent-offspring pair.
6) --headerline: Line number of the header (column names) in the vcf files.
7) --out: Path to the folder where you want to store the output file

Example:
```
python TNTallele.py --chr 22 --bim "bimpath" --maf "mafpath" --po_pairs "pairstoincludepath" --vcf "vcfpath" --headerline 4 --out "ibdmatrixpath"
```

## Step 2: Transmitted and non-transmitted allele frequencies and variance of the difference
Computing allele frequencies for the transmitted and non-transmitted alleles and the variance of their difference from the output files produced in step 1.

The following 4 arguments are required:
1) --chr: Chromosome
2) --ibdmatrix: Path to the folder where the input IBD genotype matrix (output from step 1) is stored
3) --po_pairs: Path to a file with the parent-offspring pairs you want to include in your allele frequency computations. The file should be tab delimieted, have three columns and have the header "fid ID1 ID2". Each line corresponds to one parent-offspring pair.
4) --out: Path to where you want to store the SNP wise output, IBD specific allele frequencies

Example:
```
python TNTaf1n.py --chr 22 --ibdmatrix "ibdmatrixpath" --po_pairs "pairstoincludepath" --out "ibdfreqpath"
python TNTaf1s.py --chr 22 --ibdmatrix "ibdmatrixpath" --po_pairs "pairstoincludepath" --out "ibdfreqpath"
python TNTvar.py --chr 22 --ibdmatrix "ibdmatrixpath" --po_pairs "pairstoincludepath" --out "ibdvarpath"
```

## Step 3: Make the summary statistic files
Combines the output files from step 2 in one file.
The script requires R and the following 5 arguments in this order:
1) ibdfreqpath: Path to the folder with the TNT allele frequencies (output of the previous step).
2) ibdvarpath: Path to the folder with the TNT variance (output of the previous step).
3) bim: BIM file with the SNPs that we want to include in the summary statistic file
4) maf:  Chromosome specific files containing information about MAF computed in the group of the parents for the SNPs that are in the vcf file.
5) summarystatisticpath: Path to where to store the summary statistic file

Example:
```
Rscript TNTafsumstat.R "ibdfreqpath" "ibdvarpath" "bimpath" "mafpath" "summarystatisticpath"
```

## Output file

This stepwise protocol results in a summary statistic file that has a name starting with 'TNT_AF_sumstat' and has the following columns:
1) SNP: SNP ID
2) CHR: Chromosome
3) POS: Position
4) AF: Allele frequency of the ALT allele as given in the MAF file
5) REF: Allele coded as 0
6) ALT: Allele coded as 1
7) ibd1n: Number of PO pairs for ib1nmean computations
8) ib1nmean: Allele frequency for the not-transmitted alleles
9) ibd1s: Number of PO pairs for ib1smean computations (should be the same as ibd1n).
10) ib1smean: Allele frequency for the transmitted alleles
11) ib1diffvar: Variance of the differene between the transmitted and non-transmitted allele

# Parent-offspring trios: Counting shared allele errors
Here it is shown how the numbers of errors of inferring the shared allele of double heterozygous parent-offspring pairs are counted from data of parent-offspring trios where offspring and one parent are heterozygous while the other parent is homozygous. In particular, in this count, an "error" refers to when the shared allele, inferred from phasing information, of the double heterozygous parent-offspring pair does not match with the shared allele given the homozygous genotype of the other parent. This analysis requires the ibdmatrix file from the PO analysis above as well as a file that shows the trios (unique identifiers) among the parent-offspring pairs (--po_pairs) that were used to compute the TNT allele frequencie.

## Step 1
Counting the errors. The following arguments are needed:
1) --chr: Chromosome
2) --ibdmatrix: Path to the folder where the input IBD genotype matrix (output from step 1 in the PO analysis above) is stored
3) --trios: Full path to file with PO trios (a file with columns ID1 and ID2). In this file we have one line per parent-offspring pair again. As these are trios, the identifier of each offspring will appear twice in column ID1 while each parent will appear once in the column ID2.
4) --out: Path to the folder where you want to store the output

Example:
```
python TNTtrioDH.py --chr 22 --ibdmatrix "ibdmatrixpath" --trios "triopairstoincludepath" --out "dherrorpath"
```

## Step 2
Combines the output files from step 1 in one file.
The script requires R and the following 5 arguments in this order:
1) dherrorpath: Path to the folder with the count of errors (output from step 1 above)
2) bim: BIM file with the SNPs that we want to include in the summary statistic file
3) maf:  Chromosome specific files containing information about MAF computed in the group of the parents for the SNPs that are in the vcf file.
4) summarystatisticpath: Path to where to store the summary statistic file

Example:
```
Rscript TNTtrioDHsumstat.R "dherrorpath" "bimpath" "mafpath" "summarystatisticpath"
```

## Output file

These analysis yield a summary statistic file that has a name starting with 'TNT_trio_DHerrors_' and has the following columns:
1) SNP: SNP ID
2) CHR: Chromosome
3) POS: Position
4) AF: Allele frequency of the ALT allele as given in the MAF file
5) REF: Allele coded as 0
6) ALT: Allele coded as 1
7) n_DH_other1: Number of PO trios where offspring and one parent are heterozygous for the SNP in question and the other parent is homozygous for the alternative allele (allele coded as 1).
8) n_DH_other0: Number of PO trios where offspring and one parent are heterozygous for the SNP in question and the other parent is homozygous for the reference allele (allele coded as 0).
9) error_1: Number of instances where the shared allele for the double heterozygous pairs is wrongly coded as 1 but should be 0 according to the genotype of the homozygous parent (which is 1|1).
10) error_0: Number of instances where the shared allele for the double heterozygous pairs is wrongly coded as 0 but should be 1 according to the genotype of the homozygous parent (which is 0|0).
