#!/bin/bash
cd /home/pang/Desktop/CAD_IS_PAD_20200121/data/SNP/5000
disease=$1
for i in $(seq 1 22)
do
awk -v a="${i}" '$1==a' ${disease}x >${disease}_chr${i}
awk -v a="${i}" '$1==a' txt.cytoband >cytoband${i}
done
# chr BP SNP
for t in $(seq 1 22)
do
while read p;do
# extract start
i1=$(echo $p|awk '{print $2}')
# extract end
i2=$(echo $p|awk '{print $3}')
# extract locus
i3=$(echo $p|awk '{print $6}')
# extract band
awk -v a="$i1" -v b="$i2" 'a<$2*1 && $2*1<b' ${disease}_chr${t} >${disease}_chr${t}_${i3}
# chr start end locus ann
done <cytoband${t}
done

for j in $(seq 1 22)
do
#cat chr${j}_* >chr${j}_p
awk '{print $0,FILENAME}' ${disease}_chr${j}_* >${disease}_chr${j}_p
#awk '{OFS=",";print $3,$15}' chr1_p >chr1_p0
#awk '{OFS=",";print $1,$4,$2,$3}' chr1 >chr1_l
done
cat ${disease}_chr*_p>${disease}_locus
rm ${disease}_chr*
rm cytoband*


#1mbp
for i0 in $(seq 1 22)
do
awk -v a="${i0}" '$1==a' ${disease}x >${disease}_chr${i0}
awk -v a="${i0}" '$1==a' txt.1MBP >1MBP${i0}
done
# chr BP SNP
for t0 in $(seq 1 22)
do
while read p;do
# extract start
i10=$(echo $p|awk '{print $3}')
# extract end
i20=$(echo $p|awk '{print $4}')
# extract locus
i30=$(echo $p|awk '{print $2}')
# extract band
awk -v a="$i10" -v b="$i20" 'a<$2*1 && $2*1<b' ${disease}_chr${t0} >${disease}_chr${t0}_${i30}
# chr start end locus ann
done <1MBP${t0}
done

for j0 in $(seq 1 22)
do
#cat chr${j}_* >chr${j}_p
awk '{print $0,FILENAME}' ${disease}_chr${j0}_* >${disease}_chr${j0}_p
#awk '{OFS=",";print $3,$15}' chr1_p >chr1_p0
#awk '{OFS=",";print $1,$4,$2,$3}' chr1 >chr1_l
done
cat ${disease}_chr*_p>${disease}_locus_1m
rm ${disease}_chr*
rm 1MBP*



#kb500, ld0.8 same resulty
#awk '{$4 = ($4=="False" ? $3 : $4)} 1' ld0.8_gene >ld0.8_gene_2
#awk '{print $1,$4}' ld0.8_gene_2 >ld0.8_gene
#awk '{OFS=":";print $1,$2}' ${disease}x >t
#paste ${disease}x t >${disease}c
awk 'FNR==NR{a[$1]=$0;next}($14 in a){print $0,a[$14]}' kb500_gene ${disease}c >${disease}_kb500_gene
