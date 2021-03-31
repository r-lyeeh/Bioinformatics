#!/bin/bash
# Usage: QC_sum.sh BfileName Referencefolder Workingdirectory casecontrol
file=$1
ref=$2
workdirectory=$3
# if the case and control are mixed casecontrol=1
casecontrol=$4
if [${casecontrol==1}]
then
  # Seperate case and controls
  mkdir 0_sep
  cd 0_sep
  awk '{$6==1}' ../${file}.fam >${file}_case
  awk '{$6==2}' ../${file}.fam >${file}_control
  plink --bfile ../${file} --keep ${file}_case --make-bed --out ${file}_case
  plink --bfile ../${file} --keep ${file}_control --make-bed --out ${file}_control
  ##0.1 cases
  plink --bfile ${file}_case --check-sex --out ${file}_case
  awk '{if($5=="OK") print $0}' ${file}_case.sexcheck>cases_pass_sexcheck
  plink --bfile ${file}_case --keep cases_pass_sexcheck --chr 1-22 --geno 0.05 --maf 0.01 --hwe 1e-6 --mind 0.02 --make-bed --out ${file}_case_QC
  ######flip snp to positive strand
  plink --bfile ${file}_case_QC --freq --out ${file}_case_QC
  # copy script
  cp ${workdirectory}/HRC-1000G-check-bim.pl .
  perl HRC-1000G-check-bim.pl -b ${file}_case_QC.bim -f ${file}_case_QC.frq -r ${ref}/1000GP_Phase3_combined.legend -g -p 'EUR'
  plink --bfile ${file}_case_QC --flip Strand-Flip-${file}_case_QC-1000G.txt --make-bed --out ${file}_case_QC_flip

  ##0.2 controls
  plink --bfile ${file}_control --check-sex --out ${file}_control
  awk '{if($5=="OK") print $0}' ${file}_control.sexcheck>controls_pass_sexcheck
  #wc -l controls_pass_sexcheck ###46 problems
  plink --bfile ${file}_control --keep controls_pass_sexcheck --chr 1-22 --geno 0.05 --maf 0.01 --hwe 1e-6 --mind 0.02 --make-bed --out ${file}_control_QC
  ######flip snp to positive strand
  plink --bfile ${file}_control_QC --freq --out ${file}_control_QC
  perl HRC-1000G-check-bim.pl -b ${file}_control_QC.bim -f ${file}_control_QC.frq -r ${ref}/1000GP_Phase3_combined.legend -g -p 'EUR'
  plink --bfile ${file}_control_QC --flip Strand-Flip-${file}_control_QC-1000G.txt \
  --make-bed --out ${file}_control_QC_flip

  ##################################
  ###0.3 merge cases and controls
  cd ..
  mkdir 1_SNP
  cd 1_SNP
  plink --bfile ../0_sep/${file}_case_QC_flip --bmerge ../0_sep/${file}_control_QC_flip \
  --make-bed --out ${file}_merged
  #this will produced  ${file}_merged-merge.missnp file
  #exclude these snps and re-run merging codes again
  plink --bfile ../0_sep/${file}_case_QC_flip --exclude ${file}_merged-merge.missnp --make-bed --out ${file}_case_QC_flip_1
  plink --bfile ../0_sep/${file}_control_QC_flip --exclude ${file}_merged-merge.missnp --make-bed --out ${file}_control_QC_flip_1
  ##merge using snps existed both in case and control files
  plink --bfile ${file}_case_QC_flip_1 --bmerge ${file}_control_QC_flip_1 --extract case_control_shared_snps --make-bed --out ${file}_merged
  plink --bfile ${file}_merged --missing --out ${file}_merged
else
  ###0.4 If there is only case or control
  plink --bfile ${file} --check-sex --out ${file}
  awk '{if($5=="OK") print $0}' ${file}.sexcheck>pass_sexcheck
  plink --bfile ${file} --keep pass_sexcheck --chr 1-22 --geno 0.05 --maf 0.01 --hwe 1e-6 --mind 0.02 --make-bed --out ${file}_QC
  ######flip snp to positive strand
  plink --bfile ${file}_QC --freq --out ${file}_QC
  cp ${workdirectory}/HRC-1000G-check-bim.pl .
  perl HRC-1000G-check-bim.pl -b ${file}_QC.bim -f ${file}_QC.frq -r ${ref}/1000GP_Phase3_combined.legend -g -p 'EUR'
  plink --bfile ${file}_QC --flip Strand-Flip-${file}_QC-1000G.txt --make-bed --out ${file}_merged
fi

### 1.filt SNP
###QC over the merged file
plink --bfile ${file}_merged --chr 1-22 --geno 0.05 --maf 0.01 --hwe 1e-6 --mind 0.02 --make-bed --out ${file}_merged_QC
echo "After quality check filter"
wc -l ${file}_merged_QC.fam
### 2.filt sample
cd ..
mkdir 2_sample
cd 2_sample
##IBD
plink --bfile ../1_SNP/${file}_merged_QC --genome --min 0.125 --out ${file}_merged_QC
##using sample_excluded_by_ibd.r
cp ${workdirectory}/samples_excluded_by_IBD.r .
R CMD BATCH --vanilla "--args ${workdirectory} ${file}" ./samples_excluded_by_IBD.r ./samples_excluded_by_IBD.Rout
##produce samples_excluded_by_IBD.txt
plink --bfile ../1_SNP/${file}_merged_QC --remove ${file}_samples_excluded_by_IBD.txt --make-bed --out ${file}_merged_QC_no_relative
echo "After IBD sampler filter"
wc -l ${file}_merged_QC_no_relative
###ld prunning and mds plot
mkdir population_stratification
cd population_stratification
plink --bfile ../${file}_merged_QC_no_relative --indep-pairwise 50 5 0.5 --out ${file}_merged_QC_no_relative
plink --bfile ../${file}_merged_QC_no_relative --extract ${file}_merged_QC_no_relative.prune.in --make-bed --out ${file}_merged_QC_no_relative_prune
plink --bfile ${file}_merged_QC_no_relative_prune  --cluster --mds-plot 10
cp ${workdirectory}/mds_plot.r .
R CMD BATCH --vanilla "--args ${workdirectory} ${file}" ./mds_plot.R ./mds_plot.Rout
## remove samples does not pass population stratification
cd ..
plink --bfile ${file}_merged_QC_no_relative --remove population_stratification/${file}_merged_QC_no_relative_prune.mds.sample --out ${file}_QC_all
echo "After MDS filter"
### 3.imputation
cd ..
mkdir 3_imputation
cd 3_imputation
### 3.1 remove duplicate
plink --bfile ../2_sample/${file}_QC_all --list-duplicate-vars --out ${file}
plink --bfile ../2_sample/${file}_QC_all --exclude ${file}.dupvar --make-bed --out ${file}_nodup
for i in {1..22}
do
echo "Check strand: Processing chromosome ${i}"
### 4.2 split chr
plink --bfile ${file}_nodup --chr ${i} --make-bed --out ${file}_nodup_chr${i}
### 4.3 check strand
shapeit -check -B ${file}_nodup_chr${i} -M ${ref}/genetic_map_chr${i}_combined_b37.txt -R ${ref}/1000GP_Phase3_chr${i}.hap.gz ${ref}/1000GP_Phase3_chr${i}.legend.gz ${ref}/1000GP_Phase3.sample --output-log ${file}_nodup_chr${i}.alignments
done
cp ${workdirectory}/fliplist.r
R CMD BATCH ./fliplist.r ./flipst.Rout
for i in {1..22}
do
echo "Recheck strand: Processing chromosome ${i}"
### 4.4 exclude missing
plink --bfile ${file}_nodup_chr${i} --flip flip_chr${i}.txt --make-bed --out ${file}_nodup_chr${i}_flip
### 4.5 recheck strand
shapeit -check -B ${file}_nodup_chr${i}_flip -M ${ref}/genetic_map_chr${i}_combined_b37.txt -R ${ref}/1000GP_Phase3_chr${i}.hap.gz ${ref}/1000GP_Phase3_chr${i}.legend.gz ${ref}/1000GP_Phase3.sample --output-log ${file}_nodup_chr${i}_flip.alignments
# exclude SNP
awk '{if(NR>1) print $4}' ${file}_nodup_chr${i}_flip.alignments.snp.strand> ${file}_chr${i}_exclude
plink --bfile ${file}_nodup_chr${i}_flip --exclude ${file}_chr${i}_exclude --make-bed --out ${file}_nomiss_chr${i}
done
### 4.6 exclude SNP
cp ${workdirectory}/impute_loop.sh .
for in {1..22}
do
echo "Exclude SNP: Processing chromosome ${i}"
shapeit -B ${file}_nomiss_chr${i} -M ${ref}/genetic_map_chr${i}_combined_b37.txt -R ${ref}/1000GP_Phase3_chr${i}.hap.gz ${ref}/1000GP_Phase3_chr${i}.legend.gz ${ref}/1000GP_Phase3.sample -O ${file}_chr${i}.phased -L ${file}_chr${i}.log -T 16
done
### 4.7 imputation
# Usage: impute_loop.sh start end refdir file
sbatch --nodelist=ekgen2 impute_loop.sh 1 1 ${ref} ${file}
sbatch --nodelist=ekgen3 --partition=debian impute_loop.sh 2 2 ${ref} ${file}
sbatch --nodelist=ekgen4 impute_loop.sh 3 3 ${ref} ${file}
sbatch --nodelist=ekgen5 impute_loop.sh 4 4 ${ref} ${file}
sbatch --nodelist=ekgen6 --partition=debian impute_loop.sh 5 5 ${ref} ${file}
sbatch --nodelist=ekgen7 impute_loop.sh 6 6 ${ref} ${file}
sbatch --nodelist=ekgen2 impute_loop.sh 7 7 ${ref} ${file}
sbatch --nodelist=ekgen3 --partition=debian impute_loop.sh 8 8 ${ref} ${file}
sbatch --nodelist=ekgen4 impute_loop.sh 9 9 ${ref} ${file}
sbatch --nodelist=ekgen5 impute_loop.sh 10 10 ${ref} ${file}
sbatch --nodelist=ekgen6 --partition=debian impute_loop.sh 11 12 ${ref} ${file}
sbatch --nodelist=ekgen7 impute_loop.sh 13 14 ${ref} ${file}
sbatch --nodelist=ekgen2 impute_loop.sh 15 16 ${ref} ${file}
sbatch --nodelist=ekgen3 --partition=debian impute_loop.sh 17 18 ${ref} ${file}
sbatch --nodelist=ekgen4 impute_loop.sh 19 22 ${ref} ${file}
### 4.8 merge
for i in {1..22}
do
cat ${file}_impute_chr${i}.*.out_info >${file}_merge_chr${i}_impute2.info
cat ${file}_impute_chr${i}.*.out >${file}_merge_chr${i}.impute2
sed -i "s/^---/$i/" ${file}_merge_chr${i}.impute2
### 4.9 check info
awk '{if($7>=0.9){print $2}}' ${file}_merge_chr${i}_impute2.info > ${file}_merge_chr${i}.info0.9.snps
done
### 4.10 postimpute
cp ${workdirectory}/postImputeQC.sh .
sbatch --nodelist=ekgen2 --mincpus=16 postImputQC.sh 1 1 ${file}
sbatch --nodelist=ekgen2 --mincpus=16 postImputQC.sh 2 2 ${file}
sbatch --nodelist=ekgen4 --mincpus=16 postImputQC.sh 3 3 ${file}
sbatch --nodelist=ekgen4 --mincpus=16 postImputQC.sh 4 4 ${file}
sbatch --nodelist=ekgen5 --mincpus=16 postImputQC.sh 5 5 ${file}
sbatch --nodelist=ekgen5 --mincpus=16 postImputQC.sh 6 6 ${file}
sbatch --nodelist=ekgen7 --mincpus=16 postImputQC.sh 7 7 ${file}
sbatch --nodelist=ekgen7 --mincpus=16 postImputQC.sh 8 8 ${file}
sbatch --nodelist=ekgen2 --mincpus=16 postImputQC.sh 9 9 ${file}
sbatch --nodelist=ekgen4 --mincpus=16 postImputQC.sh 10 10 ${file}
#segmentation fault for ekgen3 6 8 ??
sbatch --nodelist=ekgen5 --mincpus=16 postImputQC.sh 11 12 ${file}
sbatch --nodelist=ekgen7 --mincpus=16 postImputQC.sh 13 14 ${file}
sbatch --nodelist=ekgen2 --mincpus=16 postImputQC.sh 15 16 ${file}
sbatch --nodelist=ekgen4 --mincpus=16 postImputQC.sh 17 18 ${file}
sbatch --nodelist=ekgen5 --mincpus=16 postImputQC.sh 19 22 ${file}
