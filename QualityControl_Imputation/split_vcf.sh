#!/bin/bash
#SBATCH --mincpus=8
#SBATCH --mem=100G
file=$1
for i in {1..22};
do 
vcftools  --vcf  ${file}_QC_all.vcf  --chr $i  --recode --recode-INFO-all --out  ${file}_$i;
bgzip ${file}_${i}.recode.vcf
done
