#!/bin/bash
#Bedfiles for the individuals you want to compute pgs for
bedata=$1
#Full path of the file with the pgs weights (rsname,increasing effect allele and weight in absolute value). 
pgsfolder=$2
#Path to the folder to store the PGS
pgsout=$3
#Name of the pgs
pgs=$4

#Path to plink
plink=plink/1.90b2n/plink

$plink --bfile $bedata --score $pgsfoldergit header sum  --out $pgsout"/"$pgs
