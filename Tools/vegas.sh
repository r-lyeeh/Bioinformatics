#!/bin/bash
gene=$1
suffix=$2
folder=$3
cd ${folder}
# create snpandp
for i in $(seq 1 5)
do	
#awk '{print $2,$9}' ${gene}_${suffix}_${i}.assoc >${gene}_${suffix}_${i}.orig
vegas -G -snpandp ${gene}_${suffix}_${i}.orig -custom /home/pang/Desktop/PDE5A_Thorsten_20191011/pathway/223_4_lipids/${folder}/${gene}_${suffix}_${i} -glist /home/pang/Desktop/PDE5A_Thorsten_20191011/pathway/223_4_lipids/${gene}.glist --allow-no-sex -out ${gene}_${suffix}_${i}
awk '{print $2,$8}' ${gene}_${suffix}_${i}.out | grep -v Gene| sed 's/"//g' ${gene}_${suffix}_${i}.geneandp
vegas -P -geneandp ${gene}_${suffix}_${i}.geneandp -geneandpath /home/pang/Desktop/PDE5A_Thorsten_20191011/pathway/223_4_lipids/${gene}.vegas2pathSYM -glist /home/pang/Desktop/PDE5A_Thorsten_20191011/pathway/223_4_lipids/${gene}.glist -out ${gene}_${suffix}_${i}_path
done

