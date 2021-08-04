#!/bin/bash
#SBATCH --mincpus=16
#awk 'NR==FNR{seen1[$3]++;next};seen1[$3]==1||seen2[$3]++' t1.0 t1.0 >t2.0
for i in CAD IS PAD
do
#awk 'FNR==NR{a[$1]=$0;next}($3 in a){print $1,$2,a[$3]}' ${i}.meta snp_position_uniq >${i}_t
awk -F"\t" '{print $1}' ${i}_t >${i}_01
awk -F"\t" '{print $2,$3,$4,$5,$6,$7}' ${i}_t >${i}_02
awk '{OFS="\t";print $1,$2,$3}' ${i}_01 >${i}_03
paste -d"\t" ${i}_03 ${i}_02 >${i}_meta2
done
