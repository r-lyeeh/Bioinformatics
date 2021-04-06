#!/bin/bash
gene=$1
suffix=$2
for i in $(seq 1 5)
do
#vcftools --vcf ${gene}_${suffix}_${i}.vcf --plink --out ${gene}_${suffix}_${i}
#plink --vcf ${gene}_${suffix}_${i}.recode.vcf --recode --out ${gene}_${suffix}_${i} --double-id
plink --bfile ${gene}_${suffix}_${i} --recode --tab --out ${gene}_${suffix}_${i}
awk '{OFS=":";print "M "$1,$4}' ${gene}_${suffix}_${i}.map >${gene}_${suffix}_${i}.dat
sed -i $'1 i\\\nT CAD' ${gene}_${suffix}_${i}.dat
raremetalworker --ped ${gene}_${suffix}_${i}.ped --dat ${gene}_${suffix}_${i}.dat --traitName CAD --prefix ${gene}_${suffix}_${i}
perl calculateOddsRatio.pl ${gene}_${suffix}_${i}.CAD.singlevar.score.txt ${gene}_${suffix}_${i}.CAD.singlevar.score.OR.txt rmw
done
