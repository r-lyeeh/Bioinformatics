#!/bin/bash
working_directory=/home/pang/data/SCP_involve/GRN_20191215/Final_ld02_snpmod0
cd ${working_directory}
snplist_snp_Final=/home/pang/data/SCP_involve/GRN_20191215/h2_JohanVamsi_20190228/running_results/res_snpmod0_nonsnpcad304ld02/snplist
for chr in $(seq 1 22)
do
plinkfile_chr_ref=/raid2/Data_public/RefGenomes/1000G/release_20110521_plink/EUR.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes
plink --bfile ${plinkfile_chr_ref} --extract ${snplist_snp_Final} --make-bed --out snps_Final_chr${chr}

if [ `cat snps_Final_chr${chr}.bim | wc -l` -eq 1 ]
then
echo snps_Final_chr${chr}
awk -F"\t" '{print $2}' snps_Final_chr${chr}.bim  >> snps_allchrs_Final
fi

done
for chr in $(seq 1 22)
do
plink --bfile snps_Final_chr${chr} --indep-pairwise 500 50 0.2 --out snps_Final_chr${chr}.ld02
done
nSnp_Final_ld02=`cat snps_allchrs_Final snps_Final_chr*.ld02.prune.in | sort|uniq| wc -l`
echo ${nSnp_Final_ld02}
