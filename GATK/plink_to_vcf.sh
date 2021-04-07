#!/bin/bash
#PBS -q batch
#PBS -j oe
cd /home/pang/data/public_data/UKBB/exome_population
plink --bfile ukb_efe_chrall --recode vcf --out ukb_efe_chrall
bgzip -c ukb_efe_chrall.vcf >ukb_efe_chrall.vcf.gz
tabix -p vcf ukb_efe_chrall.vcf.gz
