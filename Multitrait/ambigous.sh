#!/bin/bash
cd /home/pang/Desktop/CAD_IS_PAD_20200121/data/SNP/ambigous_SNPi
# pick out SNPs whose OR across 1 
for i in CAD IS PAD
do 
awk '$9<1 && $10>1' ${i}_meta.meta >${i}_ambigous.meta1
awk '$9>1 && $10<1' ${i}_meta.meta >${i}_ambigous.meta2
header ${i}_meta.meta >header
cat header ${i}_ambigous.meta1 ${i}_ambigous.meta2 > ${i}_ambigous.meta
# match to meta result
awk 'FNR==NR{a[$3]=$3;next}($3 in a){print $0}' ../${i}.fdr ${i}_ambigous.meta >${i}.match
done
# IS gwas in CAD
awk 'FNR==NR{a[$3]=$3;next}($3 in a){print $0}' IS.match ../../CAD/CAD_meta.meta >ISinCAD.ambi.meta
# first OR is the correct one for IS rs2107595 is ambigous for both CAD and IS
# PAD gwas in CAD
awk 'FNR==NR{a[$3]=$3;next}($3 in a){print $0}' PAD.match ../../CAD/CAD_meta.meta >PADinCAD.ambi.meta
sort -k 2n PADinCAD.ambi.meta >t1
# first OR is the correct one for PAD;rs192231243, rs1610752, rs567589820, rs9273509, rs185913339 is PAD only
