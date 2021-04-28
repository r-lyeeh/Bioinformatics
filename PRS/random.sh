#!/bin/bash
##SBATCH --mincpus=6
##SBATCH --mem=4GB
# USAGE: replicate.sh CAD 10
trait=$1
n=$2
#workdir=/home/pang/Desktop/Liability_zeng_HS_FB_20190923/UKBB/random/bfile
for i in {1..3000}
do
#shuf -n ${n} ${workdir}/UKBB_${trait}.score >UKBB_${trait}_${n}_${i}.score
#plink --bfile ${workdir}/UKBB_${trait} --score UKBB_${trait}_${n}_${i}.score 2 3 4 sum --out UKBB_${trait}_${n}_${i}
R CMD BATCH --vanilla "--args ${trait} ${n} ${i}" random.R ./random.Rout
#awk 'BEGIN{OFS=FS="_"};NF--' ${trait}_${n}_${i}.random >t0
#mv t0 ${trait}_${n}_${i}.random
#rm t0
done
