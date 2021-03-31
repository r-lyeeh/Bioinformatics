#!/bin/bash

###################
#####post imputation QC
###################

#info>0.9
#hwe 1e-6
#maf>0.01
#genotype call reate>0.05
cd /home/pang/Desktop/ISAR_Moritz_20200910/data/1047/raw/3_imputation

chr_start=$1
chr_end=$2
file=$3
for chr in `seq $chr_start $chr_end`
do
	echo ${chr}" start"
	plink \
	--gen ${file}_merge_chr${chr}.impute2 \
	--sample ${file}_chr${chr}.phased.sample \
	--geno 0.05 --hwe 1e-6 --maf 0.01 \
	--extract ${file}_merge_${chr}.info0.9.snps \
	--make-bed \
	--out ../4_postimpute/${file}_final_chr${chr}.imputed
	echo ${chr}" end"
done
