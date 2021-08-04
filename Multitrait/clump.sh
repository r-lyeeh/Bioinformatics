#!/bin/bash
# generate clump assoc
for i in CAD IS PAD
do
awk 'FNR==NR{a[$3]=$0;next}($3 in a){print $0,a[$3]}' ../${i}meta ../${i}.fdr >${i}.fdr.meta
awk '{print $14,$15,$16,$17,$18,$21,$26}' ${i}.fdr.meta >t1
mv t1 ${i}.fdr.meta
CHR BP SNP EA NEA P BETA
