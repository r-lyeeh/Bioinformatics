#!/bin/bash
cd /home/pang/Desktop/GRN_Johan_20191215/224_mesc/UKB
source activate ldsc
#idlist=$1
GERGWAS=/home/pang/Desktop/GRN_Johan_20191215/224_mesc/UKB/ukb.assoc
#snpmodule=/home/pang/Desktop/GRN_Johan_20191215/tissue/${idlist}/snplist_module_${idlist}
#snpmodule=/home/pang/Desktop/GRN_Johan_20191215/224_mesc/snplist
for i in $(seq 1 224)
do
awk 'FNR==NR{a[$1]=$1;next}($2 in a){print $2,$4,$7}' ${snpmodule}/snplist_module_${i} ${GERGWAS} >module_${i}
sed -i '1i\\SNP A1 A2' module_${i}
#awk 'FNR==NR{a[$1]=$1;next}($2 in a){print $0}' snplist_module_all >ukbb_sumstat
#ukbsum=/home/pang/Desktop/MP_Johan_20200221/LDSC/ukbb_sumstat
# --sgined-sumstats OR,1
munge_sumstats.py --sumstats ${GERGWAS} --out module_${i}_sub --merge-alleles module_${i} --snp SNP --a1 A1 --a2 A2 --p P 

ref_ld=/home/pang/Desktop/GRN_Johan_20191215/mesc/ref_ld
w_ld=/home/pang/Desktop/GRN_Johan_20191215/mesc/ref_ld
#ldsc.py --h2 module_${i}_sub.sumstats.gz --ref-ld-chr ${ref_ld}/EUR_chr@_phase_v3_ld --w-ld-chr ${w_ld}/EUR_chr@_phase_v3_ld --out module_${i}_h2
run_mesc.py --h2med module_${i}_sub.sumstats.gz --exp-chr /home/pang/program/mesc/data/GTEx_v8/All_Tissues.@ --out module_${i}_mesc
done
grep "" module_*_mesc.all.hemed >mesc_sum.h2
sed -i -e "/Esitimate/g" mesc_sum.h2 
