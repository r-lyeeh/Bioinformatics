#!/bin/bash
#SBATCH --mincpus=16
#Usage: annotation.sh target region
trait=$1
region=$2
#sep=$3
for i in $(seq 1 22)
do
awk -v a="${i}" '$1==a' ${trait} >${trait}_chr${i}
# awk -v a="${i}" -F"," '$1==a' ${trait} >${trait}_chr${i}
awk -v a="${i}" '$1==a' ${region} >cytoband${i}
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
i3=$(echo $p|awk '{print $4}')
# extract band
awk -v a="$i1" -v b="$i2" 'a<$2*1 && $2*1<=b' ${trait}_chr${t} >${trait}_chr${t}_${i3}
#awk -v a="$i1" -v b="$i2" -F"," 'a<$2*1 && $2*1<=b' ${trait}_chr${t} >${trait}_chr${t}_${i3}
# chr start end locus ann
done <cytoband${t}
done

for j in $(seq 1 22)
do
#cat chr${j}_* >chr${j}_p
awk '{print $0,FILENAME}' ${trait}_chr${j}_* >${trait}_chr${j}_p
#awk '{OFS=",";print $3,$15}' chr1_p >chr1_p0
#awk '{OFS=",";print $1,$4,$2,$3}' chr1 >chr1_l
done
cat ${trait}_chr*_p>${trait}_${region}
rm ${trait}_chr*
rm cytoband*
awk '{print $5}' ${trait}_${region} >${trait}_t1
awk '{print $1,$2,$3,$4}' ${trait}_${region} >${trait}_t0
awk -F"_" '{OFS="_";print $3,$4}' ${trait}_t1 >${trait}_t2
#paste -d"," ${trait}_${region} ${trait}_t2 >${trait}_t3
paste -d" " ${trait}_t0 ${trait}_t2 >${trait}_${region}
rm ${trait}_t1
rm ${trait}_t2
rm ${trait}_t0
#sed -i -e 's/${trait}_chr//g' ${trait}_region

