## SNP extraction
## VCF format
cd /home/pang/Desktop/Blood_Tregouet_20200129
set=$1
filedir="/home/pang/data/public_data/PRS_basic"
vcftools --vcf ${filedir}/${set}/plink2.vcf --snps snp_id --recode --out ${set}_blood.vcf
vcftools --vcf ${filedir}/${set}/plink2.vcf --chr ${chr} --from-bp ${start} --to-bp ${end} --recode --out ${set}_pos_LDLR
#filename="$1"
#while read -r line; do
#    name="$line"
#    echo "$name"
#done < "$filename"
#extract from col 2 to col 7
awk -v b=2 -v e=7 'BEGIN{FS=OFS=" "} {for (i=b;i<=e;i++) printf "%s%s", $i, (i<e ? OFS : ORS)}' ${query}.t3
##EIDBatch
#list=`ls *.snps.recode.vcf -l | awk '{if ($5 != 0) print $9}'`
sed -i '/#.*/d' *.snps.recode.vcf
wc -l *.snps.recode.vcf >temp1.txt
sed -i '$d' temp1.txt
sed -i -e 's/.snps.recode.vcf//g' temp1.txt
awk '(NR>1) && ($1>0)' temp1.txt >temp2.txt
#rm temp1.txt
#awk 'FNR==NR{a[$2]=$1;next} ($1 in a) {print $1,$2}' temp2.txt /home/pang/data/public_data/UKBB/exome_SCP/sample.txt > EID_MAFF.txt
join -1 2 -2 1 -o 1.1,1.2,2.2 <(sort -k2 temp2.txt) <(sort -k1 /home/pang/data/public_data/UKBB/exome_SCP/sample.txt) >temp3.txt
#rm temp2.txt
awk '{a[$3]+=$1} END {for (i in a) {print i,a[i]}}' temp3.txt >MAFF_snps_sample_stat.txt
awk '{print $3}' temp3.txt |awk '{a[$1]++}END{for(i in a){print i,a[i] | "sort -r -k 2"}}' > MAFF_sample_stat.txt
#rm temp3.txt

# select random lines in a file
shuf -n N input >output
# sort -R input | head -n 100 >output

#UKBB is too large to calculate PCA or estimated grm, thus we used a pruned data set
1> Remove duplicate rs IDs
plink2 --bgen ukb_imp_chr*_QCed.bgen 'ref-first' --sample ukb_imp_chr0_QCed.sample --rm-dup 'force-first' --make-bed --out ukb_imp_chr*_QCed_rd
2> Pruned: optimizer refference from flashpca
plink2 --bfile ukb_imp_chr*_QCed_rd --indep-pairwise 1000 50 0.05 --exclude range /home/pang/program/flashpca/exclusion_regions_hg19.txt --out ukb_imp_chr*_QCed_rd
plink2 --bfile ukb_imp_chr*_QCed_rd --extract ukb_imp_chr*_QCed_rd.prune.in --make-bed --out ukb_imp_chr*
3> GRM extimated
gcta64 --bfile ukb_imp_chr* --make-grm-bin --out ukb_imp_chr*
4> Calculate PCA
gcta64 --mgrm multi_grm_ukb_imp.txt --make-grm-bin --out ukb_imp
gcta64 --grm ukb_imp --pca 20 --thread-num 4 --out ukb_imp_20PSs
