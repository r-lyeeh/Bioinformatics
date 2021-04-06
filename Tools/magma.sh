#!/bin/bash
cd /home/pang/Desktop/MAFF_Moritz_20190927/PRS/PRS_OR/test/G7/magma
#convert vcf to bfile
#set=$1
#gene=$2
#cd ${set}
#plink2 --vcf ${set}_pos_${gene}.recode.vcf --make-bed --out ${set}_${gene}
#magma
gene=$1
geneLoci="/home/pang/program/magma/ref/NCBI37.3.gene.loc" 
#magma --annotate --snp-loc ${set}_${gene}.bim --gene-loc ${geneLoci} --out ${set}_${gene}.genes.annot
#magma --bfile ${set}_${gene} --gene-annot ${set}_${gene}.genes.annot --out ${set}_${gene}_magma
for i in $(seq 1 5) 
do
#plink2 --vcf ${gene}_LDL_${i}.recode.vcf --make-bed --out ${gene}_LDL_${i}
magma --annotate --snp-loc ${gene}_LDL_1.bim --gene-loc ${geneLoci} --out ${gene}_LDL_${i}
magma --bfile ${gene}_LDL_${i} --gene-annot ${gene}_LDL_${i}.genes.annot --out ${gene}_LDL_${i}_magma
done
