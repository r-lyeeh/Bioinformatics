#!/bin/bash
#ln -d /home/pang/Desktop/GWAS_all/scripts/mergelist
file=$1
code=$2
## unzip michgan server result
for i in {1..22}
do
#IR5 S0OdPy0hT<vi+B
#case_munich bEK1c3qCCujRV
#Procam2 1HgPo7ijTqwZQP
#cotrol_focus 'woL9ZnoAPMwTu4'
unzip -P ${code} chr_${i}.zip
gunzip chr${i}.info.gz 
plink --const-fid --vcf chr${i}.dose.vcf.gz --make-bed --out chr${i}.dose 
## postimpute
awk '{if($6>=0.9){print $1}}' chr${i}.info >chr${i}.info0.9.snps
plink --bfile chr${i}.dose --geno 0.05 --hwe 1e-6 --maf 0.01 --extract chr${i}.info0.9.snps --make-bed --out chr${i}.post
done
## convert vcf to bfile
plink --bfile chr1.post --merge-list mergelist --make-bed --out ${file}_HRC
