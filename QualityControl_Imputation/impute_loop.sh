#!/bin/bash
chr_from=$1
chr_to=$2
refdir=$3
file=$4
for ((j=${chr_from}; j<=${chr_to}; j++))
	do
	for i in {0..49}
		do
		c=$(($i*5))
		start="${c}000001"
		end="$(($c+5))000000"
		impute2 -h ${refdir}/1000GP_Phase3_chr${j}.hap.gz \
		-l ${refdir}/1000GP_Phase3_chr${j}.legend.gz \
		-m ${refdir}/genetic_map_chr${j}_combined_b37.txt \
		-use_prephased_g \
		-known_haps_g ${file}_${j}.phased.haps \
		-int $start $end -buffer 500 -pgs_miss -filt_rules_l 'EUR==0' \
		-o ${file}_impute_chr${j}.$(($i+1)).out
		done
	done
