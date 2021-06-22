#!/bin/bash
#SBATCH --mincpus=12
#SBATCH --mem=8GB
# USAGE: replicate.sh CAD 10
trait=$1
n=$2
for i in {1..3000}
do
#shuf -n ${n} UKBB_${trait}.score >UKBB_${trait}_${n}_${i}.score
#plink --bfile UKBB_${trait} --score UKBB_${trait}_${n}_${i}.score 2 3 4 sum --out UKBB_${trait}_${n}_${i}
#using plink to calculate the PRS, needs to specify which is the ture risk allele first
R CMD BATCH --vanilla "--args ${trait} ${n} ${i}" smote_random.R ./smote_random.Rout
#awk 'BEGIN{OFS=FS="_"};NF--' ${trait}_${n}_${i}.random >t0
#mv t0 ${trait}_${n}_${i}.random
#rm t0
done
