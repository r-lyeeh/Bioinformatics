#!/bin/bash
## Usage: bedmap.sh SNP Gene range outname
## SNP/Gene format: chr start end SNP/Gene
SNP=$1
Gene=$2
range=$3
outname=$4
for i in $(seq 1 22)
do
awk -v a="$i" '$1==a' ${SNP} >t1.2
awk -v a="$i" '$1==a' ${Gene} >t1.1
bedmap --echo --echo-map-id-uniq --range ${range} t1.2 t1.1 >res_${i}
#closest-features --closest --dist t1.1 t1.2 >res_${i}
done
cat res_* >${outname}
#cat res_* >gene_100kb_fourier
rm t1.*
rm res_*
