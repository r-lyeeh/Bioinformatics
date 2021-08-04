#!/bin/bash
#SBATCH --mem=70G
#SBATCH --mincpus=16
cd /home/pang/Desktop/CAD_IS_PAD_20200121/data/METAL/gene_annotation

molocres=$1
awk '{print $2,$4}' ${molocres} >t00
awk 'FNR==NR{a[$2]=$0;next}($3 in a){print $0,a[$3]}' t00 /home/pang/Desktop/CAD_IS_PAD_20200121/data/METAL/gene_annotation/snp_position_uniq_tab >t0
awk '{print $1,$2,$2+1,$3,$4}' t0 >t0.1
sort -n -k1,1 -k2,2 t0.1 >t1
for i in $(seq 1 22)
do
awk -v a="$i" '$1==a' t1 >t1.1
#sort -n -k 2 t1 >t1.1
awk -v a="$i" '$1==a' NCBI.gene >t1.2
bedmap --echo --echo-map-id-uniq --range 500000 t1.1 t1.2 >res_${i}
#closest-features --closest --dist t1.1 t1.2 >res_${i}
done
cat res_* >gene_500kb_500kb
#cat res_* >gene_100kb_fourier
#rm t0*
#rm t1*
