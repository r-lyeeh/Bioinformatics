#!/bin/bash
work_directory=$1
cd ${work_directory}
snp=$2
# Region annotation
for i in $(seq 1 22)
do
#awk -v a="$i" '$1==a' loci_region_unique >chr${i}
awk -v a="chr${i}" '$1==a' snp_${snp} >chr${i}
awk -v a="chr${i}" '$1==a' cytoBand.txt >cytoband${i}
#awk -v a="chr${i}" '$1==a' 1MBP.txt >cytoband${i}
done

for t in $(seq 1 22)
do
while read p;do
# extract chr
i1=$(echo $p|awk '{print $1}')
# extract start
i2=$(echo $p|awk '{print $2}')
# extract end
i3=$(echo $p|awk '{print $3}')
# extract band
i4=$(echo $p|awk '{print $4}')
awk -v a="$i1" -v b="$i2" 'a<$2*1 && $2*1<b' chr${t} >chr${t}_${i3}
#awk -v a="$i1" -v b="i2" 'a<$2*1 && $2*1<b' ukb_chr${t}.phenotype.glm.logistic.hybrid >chr${t}_${i4}
# must *1!!!!! Otherwise the type is not numeric
done <chr${t}
done

for j in $(seq 1 22)
do
#cat chr${j}_* >chr${j}_p
awk '{print $0,FILENAME}' chr${j}_* >chr${j}_p
#awk '{OFS=",";print $3,$15}' chr1_p >chr1_p0
#awk '{OFS=",";print $1,$4,$2,$3}' chr1 >chr1_l
done
cat chr*_p>locus_1m_${snp}
#cat chr*_p>locus_${snp}
rm chr*
rm cytoband*
awk '{print $4}' locus_1m_${snp} |sort |uniq >t1
