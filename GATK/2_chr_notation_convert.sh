#!/bin/bash
#PBS -q batch
#PBS -j oe
for i in {1..22} 
do
echo "$i chr$i" >> chr_name_conv.txt
done
echo "23 chrX" >>chr_name_conv.txt
echo "24 chrY" >>chr_name_conv.txt
bcftools annotate --rename-chrs chr_name_conv.txt ukb_efe_chrall.HC.vcf2.gz -Oz -o ukb_efe_chrall.HC.vcf3.gz
plink --bfile ukb_efe_chrall --recode vcf --out ukb_efe_chrall
bgzip -c ukb_efe_chrall.vcf >ukb_efe_chrall.vcf.gz
tabix -p vcf ukb_efe_chrall.vcf.gz

