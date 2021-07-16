#!/bin/bash
workingdirectory=pwd
# raw QC
plink --bfile ${file} --chr 1-22 --geno 0.05 --maf 0.01 --hwe 1e-6 --mind 0.02 --make-bed --out ${file}_QC
# filt sample
plink --bfile ${file}_QC --genome --min 0.125 --out ${file}_merged_QC
##using sample_excluded_by_ibd.r
R CMD BATCH --vanilla "--args ${workdirectory} ${file}" ./samples_excluded_by_IBD.r ./samples_excluded_by_IBD.Rout
##produce samples_excluded_by_IBD.txt
plink --bfile ${file}_QC --remove ${file}_samples_excluded_by_IBD.txt --make-bed --out ${file}_merged_QC_no_relative
###ld prunning and mds plot
mkdir population_stratification
cd population_stratification
plink --bfile ../${file}_merged_QC_no_relative --indep-pairwise 50 5 0.5 --out ${file}_merged_QC_no_relative
plink --bfile ../${file}_merged_QC_no_relative --extract ${file}_merged_QC_no_relative.prune.in --make-bed --out ${file}_merged_QC_no_relative_prune
plink --bfile ${file}_merged_QC_no_relative_prune  --cluster --mds-plot 10
R CMD BATCH --vanilla "--args ${workdirectory} ${file}" ./mds_plot.R ./mds_plot.Rout
## remove samples does not pass population stratification
cd ..
plink --bfile ${file}_merged_QC_no_relative --remove population_stratification/${file}_merged_QC_no_relative_prune.mds.sample --make-bed --out ${file}_QC_all
######flip snp to positive strand
plink --bfile ${file}_QC_all --freq --out ${file}_QC_all
perl HRC-1000G-check-bim.pl -b ${file}_QC_all.bim -f ${file}_QC_all.frq -r ${ref}/1000GP_Phase3_combined.legend -g -p 'EUR'
perl HRC-1000G-check-bim.pl -b ${file}_QC_all.bim -f ${file}_QC_all.frq -r ~/data/public_data/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -p 'EUR'
plink --bfile ${file}_QC_all --flip Strand-Flip-${file}_control_QC-1000G.txt --make-bed --out ${file}_QC_all_flip
## create vcf
plink --bfile ${file}_QC_all_flip --recode vcf --out ${file}_QC_all
## seperate vcf
sh split_vcf.sh
## upload to michgan server
python mychgan.py
