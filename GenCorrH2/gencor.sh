#!/bin/bash
sumstat1=$1
sumstat2=$2
LDscores=$3
rg=$4

ldsc.py \
--rg $sumstat1,$sumstat2 \
--ref-ld $LDscores \
--w-ld $LDscores \
--intercept-gencov 0,0 \
--out $rg"_rg"
