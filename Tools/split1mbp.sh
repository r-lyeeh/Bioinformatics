#!/bin/bash
cd /home/pang/data/public_data/SNP/locus
snp=$1
for i in $(seq 1 22)
do
awk -v a="chr${i}" '$1==a' snp_${snp} >chr${i}
awk -v a="chr${i}" '$1==a' 1MBP.txt >cytoband${i}
done
# chr BP SNP
for t in $(seq 1 22)
do
while read p;do
# extract start
i1=$(echo $p|awk '{print $3}')
# extract end
i2=$(echo $p|awk '{print $4}')
# extract locus
i3=$(echo $p|awk '{print $2}')
# extract band
awk -v a="$i1" -v b="$i2" 'a<$2*1 && $2*1<b' chr${t} >chr${t}_${i3}
# chr start end locus ann
done <cytoband${t}
done

for j in $(seq 1 22)
do
#cat chr${j}_* >chr${j}_p
awk '{print $0,FILENAME}' chr${j}_* >chr${j}_p
#awk '{OFS=",";print $3,$15}' chr1_p >chr1_p0
#awk '{OFS=",";print $1,$4,$2,$3}' chr1 >chr1_l
done
cat chr*_p>locus_1m_${snp}
rm chr*
rm cytoband*
awk '{print $4}' locus_1m_${snp} |sort |uniq >t1
