#!/bin/bash
#SBATCH --nodelist=ekgen5
#SBATCH --job-name=tissue_full_CT
#SBATCH --output=tissue_full_CT-%j.out

### get hsq for snplist full list
working_directory=/home/pang/Desktop/GRN_Johan_20191215/MP
cd ${working_directory}

name_out=snpmod_full_20200102

echo name_snplist Variance_VgVgL Se_VgVgL Pval \
nSnp_Orig nSnp_OrigInRef nSnp_OrigInRef_ld05 nSnp_OrigInRef_ld02 \
nSnp_Avail nSnp_Avail_ld05 nSnp_Avail_ld02 \
nSnp_Final nSnp_Final_ld05 nSnp_Final_ld02 > hsq_res_${name_out}
oldIFS="$IFS"
IFS=$'\n' arr=($(</home/pang/Desktop/GRN_Johan_20191215/MP/module))
IFS="$oldIFS"
b=${#arr[@]}
c=$((b-1))
echo "${b}"
for i in $(seq 0 ${c}) 
do
id_snplist=${arr[${i}]}
orig_snplist=${working_directory}/${id_snplist}/snplist_module_${id_snplist}
name_snplist=snpmod${id_snplist}
echo ${name_snplist}

snplist=${orig_snplist}
snplist_snp_Avail=${working_directory}/${id_snplist}/res_${name_snplist}_fullsnp/snp_Avail
snplist_snp_Final=${working_directory}/${id_snplist}/res_${name_snplist}_fullsnp/snplist

# hsq result:
hsqfile=${working_directory}/${id_snplist}/res_${name_snplist}_fullsnp/gcta/all9study_Imput_snpmod${id_snplist}.hsq
Variance_VgVgL=`cat ${hsqfile} | grep 'V(G)/Vp_L' | awk -F"\t" '{print $2}'`
Se_VgVgL=`cat ${hsqfile} | grep 'V(G)/Vp_L' | awk -F"\t" '{print $3}'`
Pval=`cat ${hsqfile} | grep 'Pval' | awk -F"\t" '{print $2}'`

echo  ${Variance_VgVgL} ${Se_VgVgL} ${Pval}


# not pruned
nSnp_Orig=`cat ${snplist} | wc -l`
nSnp_Avail=`cat ${snplist_snp_Avail} | wc -l`
nSnp_Final=`cat ${snplist_snp_Final} | wc -l`
echo ${nSnp_Orig} ${nSnp_Avail} ${nSnp_Final}


# nSnp ld-independent
mkdir temp_${name_out}
cd temp_${name_out}

# - list of all snps
> snps_allchrs_AvailInRef
> snps_allchrs_Avail
> snps_allchrs_Final
for chr in $(seq 1 22)
do
plinkfile_chr_ref=/raid2/Data_public/RefGenomes/1000G/release_20110521_plink/EUR.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes

# available in Reference
plink --bfile ${plinkfile_chr_ref} --extract ${snplist} --make-bed --out snps_AvailInRef_chr${chr}

if [ `cat snps_AvailInRef_chr${chr}.bim | wc -l` -eq 1 ]
then
echo snps_AvailInRef_chr${chr}
awk -F"\t" '{print $2}' snps_AvailInRef_chr${chr}.bim  >> snps_allchrs_AvailInRef
fi

# available in analysis
plink --bfile ${plinkfile_chr_ref} --extract ${snplist_snp_Avail} --make-bed --out snps_Avail_chr${chr}

if [ `cat snps_Avail_chr${chr}.bim | wc -l` -eq 1 ]
then
echo snps_Avail_chr${chr}
awk -F"\t" '{print $2}' snps_Avail_chr${chr}.bim  >> snps_allchrs_Avail
fi

# final snps in analysis
plink --bfile ${plinkfile_chr_ref} --extract ${snplist_snp_Final} --make-bed --out snps_Final_chr${chr}

if [ `cat snps_Final_chr${chr}.bim | wc -l` -eq 1 ]
then
echo snps_Final_chr${chr}
awk -F"\t" '{print $2}' snps_Final_chr${chr}.bim  >> snps_allchrs_Final
fi

done



# - list of pruned snps
# do pruning chr by chr
for chr in $(seq 1 22)
do

# ld05prune
plink --bfile ./snps_AvailInRef_chr${chr} --indep-pairwise 500 50 0.5 --out snps_AvailInRef_chr${chr}.ld05

plink --bfile ./snps_Avail_chr${chr} --indep-pairwise 500 50 0.5 --out snps_Avail_chr${chr}.ld05

plink --bfile ./snps_Final_chr${chr} --indep-pairwise 500 50 0.5 --out snps_Final_chr${chr}.ld05

# ld02prune
plink --bfile ./snps_AvailInRef_chr${chr} --indep-pairwise 500 50 0.2 --out snps_AvailInRef_chr${chr}.ld02

plink --bfile ./snps_Avail_chr${chr} --indep-pairwise 500 50 0.2 --out snps_Avail_chr${chr}.ld02

plink --bfile ./snps_Final_chr${chr} --indep-pairwise 500 50 0.2 --out snps_Final_chr${chr}.ld02

done


# - get results of nsnps
nSnp_OrigInRef=`cat ./snps_AvailInRef_chr*.bim | wc -l`
echo  ${nSnp_OrigInRef}

nSnp_OrigInRef_ld05=`cat ./snps_allchrs_AvailInRef ./snps_AvailInRef_chr*.ld05.prune.in | sort|uniq | wc -l`
nSnp_Avail_ld05=`cat ./snps_allchrs_Avail ./snps_Avail_chr*.ld05.prune.in | sort|uniq| wc -l`
nSnp_Final_ld05=`cat ./snps_allchrs_Final ./snps_Final_chr*.ld05.prune.in | sort|uniq| wc -l`
echo ${nSnp_OrigInRef_ld05} ${nSnp_Avail_ld05} ${nSnp_Final_ld05}


nSnp_OrigInRef_ld02=`cat ./snps_allchrs_AvailInRef ./snps_AvailInRef_chr*.ld02.prune.in | sort|uniq | wc -l`
nSnp_Avail_ld02=`cat ./snps_allchrs_Avail ./snps_Avail_chr*.ld02.prune.in | sort|uniq| wc -l`
nSnp_Final_ld02=`cat ./snps_allchrs_Final ./snps_Final_chr*.ld02.prune.in | sort|uniq| wc -l`
echo ${nSnp_OrigInRef_ld02} ${nSnp_Avail_ld02} ${nSnp_Final_ld02}


cd ..
# rm -r temp_${name_out}


echo ${name_snplist} ${Variance_VgVgL} ${Se_VgVgL} ${Pval} \
${nSnp_Orig} ${nSnp_OrigInRef} ${nSnp_OrigInRef_ld05} ${nSnp_OrigInRef_ld02} \
${nSnp_Avail} ${nSnp_Avail_ld05} ${nSnp_Avail_ld02} \
${nSnp_Final} ${nSnp_Final_ld05} ${nSnp_Final_ld02}  >> hsq_res_${name_out}
done


cat hsq_res_${name_out}
