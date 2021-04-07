#!/bin/bash
#SBATCH --nodelist=ekgen7
#SBATCH --job-name=ldsc
#SBATCH --output=ldsc-%j.out

cd /home/pang/Desktop/MP_Johan_20200221/LDSC
oldIFS="$IFS"
IFS=$'\n' arr=($(</home/pang/Desktop/MP_Johan_20200221/snp_module/module))
IFS="$oldIFS"
b=${#arr[@]}
c=$((b-1))
echo "${b}"
for i in $(seq 0 ${c}) 
do
id_snplist=${arr[${i}]}
sbatch --nodelist=ekgen7 --job-name=ldsc_${id_snplist} --output=ldsc_${id_snplist}-%j.out ldsc_excutes.sh ${id_snplist}
done

