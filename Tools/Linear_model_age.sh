#!/bin/bash
cd /home/pang/Desktop/Scattered/NLRP3/vcf2
#for i in $(seq 2 7)
#do
#vcftools --vcf G${i}.vcf --snp rs10754555 --recode --out G${i}_NLRP3
#done
awk '$6>=60' merge_rs12048215.clean_1.fam >merge_NLRP3_over_60_v1.fam
awk '$6<60' merge_rs12048215.clean_1.fam >merge_NLRP3_under_60_v1.fam
awk '$6>=50' merge_rs12048215.clean_1.fam >merge_NLRP3_over_50_v1.fam
awk '$6<50' merge_rs12048215.clean_1.fam >merge_NLRP3_under_50_v1.fam
for i in under_60 over_60 under_50 over_50
do
# crude
plink --bfile merge_rs12048215.clean --keep merge_NLRP3_${i}_v1.fam --out merge_NLRP3_${i} --make-bed
awk '{print $1,$2,$3,$4,$5,$7}' merge_NLRP3_${i}_v1.fam >merge_NLRP3_${i}.fam
plink --bfile merge_NLRP3_${i} --ci 0.95 --assoc --out merge_NLRP3_${i}
sed -i '1i\\FID IID F1 F2 SEX AGE CAD' merge_NLRP3_${i}_v1.fam
#genotypic
plink --bfile merge_NLRP3_${i} --ci 0.95 --logistic genotypic --out merge_NLRP3_${i}_add_una
#adjusted_sex_age_genotypic
plink --bfile merge_NLRP3_${i} --ci 0.95 --pheno merge_NLRP3_${i}_v1.fam --pheno-name CAD --covar merge_NLRP3_${i}_v1.fam --covar-name AGE,SEX --logistic genotypic --out merge_NLRP3_${i}_add_adj
#adjusted_sex_age
plink --bfile merge_NLRP3_${i} --ci 0.95 --pheno merge_NLRP3_${i}_v1.fam --pheno-name CAD --covar merge_NLRP3_${i}_v1.fam --covar-name AGE,SEX --logistic --out merge_NLRP3_${i}_adj
#sum_score
plink --bfile merge_NLRP3_${i} --recode vcf --out merge_NLRP3_${i}
bgzip merge_NLRP3_${i}.vcf
bcftools stats -s- merge_NLRP3_${i}.vcf.gz >merge_NLRP3_${i}_bcftools.stats
# extract sum score
grep -hnr "PSC" merge_NLRP3_${i}_bcftools.stats >merge_NLRP3_${i}.sumscore
awk '{print $3,$4,$5,$6,$14}' merge_NLRP3_${i}.sumscore >temp1
sed -i '/sample/d' temp1
sed -i '1isample nRefHom nNonRefHom nHets nMissing' temp1
mv temp1 merge_NLRP3_${i}.sumscore
done
