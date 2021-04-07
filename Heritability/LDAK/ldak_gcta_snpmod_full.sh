#!/bin/bash
#SBATCH --mincpus=16
#SBATCH --partition=debian
#SBATCH --export=ALL
id_snplist=$1

n_thd=6


working_directory=/home/pang/Desktop/GRN_Johan_20191215
cd ${working_directory}


# - snplist
snplist=${working_directory}/${id_snplist}/snplist_module_${id_snplist}
name_snplist=snpmod${id_snplist}
echo ${snplist}
echo ${name_snplist}


# - individual-level genotypes
folder_indsgeno=/raid2/zeng/work_zeng_summary/work_zeng_as_main_author/20150624-20190221_h2Johan/Heritability_Johan_20150624/prep_20151023/pre_h2_9study_20151023/


# run fullsnp here! i.e. without exclusion of snp_gwascad
### input & run ~~~~~~~~~~~~~~~~~~

mkdir ${working_directory}/${id_snplist}/res_${name_snplist}_fullsnp
cd ${working_directory}/${id_snplist}/res_${name_snplist}_fullsnp

this_snplist=${snplist}

# availability of original snplist
# - based on individual-level genotypes
name_prefix=all9study_Imput

for chr in $(seq 1 22)
do
	plinkfile_chr=${folder_indsgeno}plinkfile_${name_prefix}/${name_prefix}_chr${chr}
	plink --bfile ${plinkfile_chr} --extract ${this_snplist} --make-bed --out availin_chr${chr}_${name_snplist}
done
rm availin_chr*_${name_snplist}.log
rm availin_chr*_${name_snplist}.bed
rm availin_chr*_${name_snplist}.fam
# only the available SNP-IDs are needed (.bed, .fam not needed)
# how many SNPs available in the all post-imputation-QCed data, here show # SNPs
cat availin_chr*_${name_snplist}.bim | wc -l
# how this snplist span across 22 chromosomes, here show # chrs
ls availin_chr*_${name_snplist}.bim | cut -d'_' -f2  | wc -l
ls availin_chr*_${name_snplist}.bim | cut -d'_' -f2  > chr_Avail
# prepare a list of SNPs not found in the prepared post-imputation-QCed data
cat availin_chr*_${name_snplist}.bim | awk -F"\t" '{print $2}' > snp_Avail
cat ${this_snplist} ./snp_Avail | sort | uniq -u > snp_NotAvail

### 0 prepare data
current_path=`pwd`
for chr in $(seq 1 22)
do
	# find proxy-SNP
	plinkfile_chr_ref=/raid2/Data_public/RefGenomes/1000G/release_20110521_plink/EUR.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes

	# get an extended set of SNPs based on reference by search for the up/down 500kb flankings and with ld > 0.8
	plink --bfile ${plinkfile_chr_ref} --ld-snp-list snp_NotAvail --r2 --ld-window-kb 500 --ld-window 19999 --ld-window-r2 0.8 --out snp_NotAvail_chr${chr}_500kbld08

rm snp_NotAvail_chr${chr}_500kbld08.nosex
rm snp_NotAvail_chr${chr}_500kbld08.log

R CMD BATCH --vanilla "--args snp_NotAvail_chr${chr}_500kbld08.ld ${current_path}" ${working_directory}/find_proxy_snp.R

done

cat snp_NotAvail_BestProxy_chr* | awk -F" " '{print $6}'  > snp_NotAvail_BestProxyFound
ls snp_NotAvail_BestProxy_chr* |cut -d'_' -f4 > chr_NotAvail_BestProxyFound


cat snp_NotAvail_chr*_500kbld08.ld | awk -F" " '{print $6}' | sort > snp_NotAvail_500kbld08
sed -i '/^SNP_B/d'  snp_NotAvail_500kbld08
ls snp_NotAvail_chr*_500kbld08.ld |cut -d'_' -f3 > chr_NotAvail_500kbld08


###
mkdir prepdata
mv availin_* prepdata/
mv *.ld prepdata/
mv find_proxy_snp.Rout prepdata/
mv snp_NotAvail_BestProxy_chr* prepdata/




### snp for analysis:
cat snp_Avail snp_NotAvail_BestProxyFound | sort | uniq > snplist
cat chr_Avail chr_NotAvail_BestProxyFound | sort | uniq | sed 's/^chr//g' > chrlist

#cat snp_Avail snp_NotAvail_500kbld08 | sort | uniq > snplist
#cat chr_Avail chr_NotAvail_500kbld08 | sort | uniq | sed 's/^chr//g' > chrlist

n_chr=`cat ./chrlist | wc -l`


###
mkdir ./plinkfile
cd  ./plinkfile

for i_chr in $(seq 1 ${n_chr})
do
    chr=`sed -n ${i_chr}p ../chrlist`

	# prepare for analysis based on individual-level genotypes
	plinkfile_chr=${folder_indsgeno}plinkfile_${name_prefix}/${name_prefix}_chr${chr}
	plink --bfile ${plinkfile_chr} --extract ../snplist --make-bed --out chr${chr}_${name_snplist}

	# also prepare sample id files for cases and controls
	awk -F" " '{if($6==1) print $1" "$2}' ${plinkfile_chr}.fam > sub1
	awk -F" " '{if($6==2) print $1" "$2}' ${plinkfile_chr}.fam > sub2

done
#rm chr*_${name_snplist}.log
cat chr*_${name_snplist}.bim | wc -l

#memo original list:
cat ${snplist} | wc -l


### 1. calculate SNP weights using LDAK (adjusting for LD)
mkdir ../ldak
cd ../ldak


for i_chr in $(seq 1 ${n_chr})
do
	chr=`sed -n ${i_chr}p ../chrlist`
	this_plinkfile=../plinkfile/chr${chr}_${name_snplist}
# 1. calculate weights
# 	it is no longer possible to calculate weights twice; a better solution is to first thin predictors to remove duplicates
#	(which LDAK now does by default when cutting predictors into sections)
	ldak5.linux.fast --cut-weights sections_chr${chr} --bfile ${this_plinkfile}
	# now The details of each section will be stored in sections_chr/section.details

	#Compute weightings for each SNP section
	#The exact number of sections is provided in sections/section.number
	nl=`cat sections_chr${chr}/section.details | wc -l`
	n_sct=`expr ${nl} - 4`
	echo n_sct = ${n_sct}

	# calculate weights for cases and controls separately (called 'Subset' in ldak5) , and then join.
	for i in $(seq 1 ${n_sct})
	do
	ldak5.linux.fast --calc-weights sections_chr${chr}  --bfile ${this_plinkfile} --subset-number 2 --subset-prefix ../plinkfile/sub --section ${i}
	echo calc-weights section_1st ${i} finished.
	done
	ldak5.linux.fast --join-weights sections_chr${chr} --bfile ${this_plinkfile}
	# The final weightings will be stored in sections_chr${chr}/weights.all.

# 2. Compute kinships for this chromosome
	# 	get a list of SNP position
	#	awk '{print $4 > "chr"$1""}' ${this_plinkfile}.bim
	#	ldak5.beta.fast --cut-kins kinships_chr${chr} --bfile ${this_plinkfile} --by-chr YES
	#	ldak5.beta.fast --calc-kins kinships_chr${chr} --bfile ${this_plinkfile} --partition 1 --weights sections_chr${chr}/weights.all --power -0.25
	#	you must use "--power" to specify how to standardize predictors; the old default was -1, but we now suggest -0.25
	#	ldak5.beta.fast --calc-kins-direct ./chr${chr}_${name_snplist} --bfile ${this_plinkfile} --weights sections_chr${chr}/weights.all --power -0.25
	# ! :( unfortunately, the ldak5 won't generate grm.N.bin!  , which is nachher essential for gcta!, so use ldak4 instead here

	nweightsfile=`cat sections_chr${chr}/weights.all | wc -l`
	nsnp=`expr ${nweightsfile} - 1`
	awk -F" " -v nsnp=${nsnp} '{if(NR>1) print $2" "nsnp" "nsnp" "nsnp" "$1}' sections_chr${chr}/weights.all > sections_chr${chr}/weightsALL
	ldak.out --calc-kins-direct ./chr${chr}_${name_snplist} --bfile ${this_plinkfile} --weights sections_chr${chr}/weightsALL
done

#
echo "chr-wise GRM calculation done by ldak5 + ldak4 ."


###
mkdir ../gcta
cd ../gcta


# 2. merge multiple grm

# snp for analysis:
# ../chrlist
# ../snplist
# however, final plinkfile only contain SNPs available
ls ../plinkfile/chr*_${name_snplist}.bim | cut -d'_' -f1  > a
cat a |cut -d'/' -f3 | sort | uniq | sed 's/^chr//g'  > chrlist_final
n_chr=`cat ./chrlist_final | wc -l`

> multi_grm_${name_prefix}_${name_snplist}.txt
for i_chr in $(seq 1 ${n_chr})
do
chr=`sed -n ${i_chr}p chrlist_final`
echo ../ldak/chr${chr}_${name_snplist} >> multi_grm_${name_prefix}_${name_snplist}.txt
done

gcta64  --mgrm multi_grm_${name_prefix}_${name_snplist}.txt  --make-grm-bin --out ${name_prefix}_${name_snplist} --thread-num ${n_thd}




# 3. calculate h2,  gcta

# - gcta formated phenotype
# 3columns: FID IID plinkpheno (case-2, control-1.)
pheno_gcta=${folder_indsgeno}all9study.phen

# - gcta formated PCs
qcovar_gcta=${folder_indsgeno}all9study_Geno_20PCs.txt.eigenvec



# now calculate h2 via reml
gcta64  --reml --grm ${name_prefix}_${name_snplist} --prevalence 0.05 --pheno ${pheno_gcta} --qcovar ${qcovar_gcta} --out ${name_prefix}_${name_snplist} --thread-num ${n_thd}
