#!/bin/bash
#Full path to summary statistic file
insumstat=$1
#Folder path and prefix to LD scores
LDscores=$2
#Folder path and prefix to the output summary statistic file
outpath=$3

#Heritability constraining intercept to 1
ldsc.py \
--h2 $insumstat \
--ref-ld $LDscores \
--intercept-h2 1 \
--w-ld $LDscores \
--out $outpath"_h2"
