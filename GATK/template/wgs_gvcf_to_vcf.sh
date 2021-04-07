#!/bin/bash
#fastq_to_gvcf

#path
gatk=/home/pang/program/gatk/gatk

#reference
reference=/home/pang/program/gatk/GATK_bundle/Homo_sapiens_assembly38.fasta
GATK_bundle=/home/pang/program/gatk/GATK_bundle

## shell parameter
samples=$1 #sample ID seperated by ","
#file=/home/pang/data/samples,txt
#while IFS= read -r line; do
#	sample_gvcfs=${sample_gvcfs}"-V $indir/${line}_23161_0_0.gvcf.gz "
#done <"$file"
indir=$2 #same as the outdir of fastq to gvcf
outname=$3 #prefix of out file

outdir=$indir

#population variation outdir
if [ ! -d $outdir/population ]
then mkdir -p $outdir/population
fi

#extract sample ID
samples=$(echo $samples | tr "," "\n")

#Joint genotyping
#joint - genotypeGVCFG
sample_gvcfs=""
for sample in $samples ; do
	sample_gvcfs=${sample_gvcfs}"-V $outdir/${sample}/gatk/${sample}.HC.g.vcf.gz "
done
time $gatk CombineGVCFs \
	-R $reference/Homo_sapiens_assembly38.fasta \
	${sample_gvcgs} \
	-O $outdir/population/${outname}.HC.g.vcf.gz && echo "** ${outname}.HC.g.vcf.gz done **"
time $gatk GenotypeGVCFs \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-V $outdir/population/${outname}.HC.g.vcf.gz \
	-O $outdir/population/${outname}.HC.vcf.gz && echo "** ${outname}.HC.vcf.gz done **"
#joint chr - genotype
#chrom=(chr1 chr2 chr3 chr4...
#for i in ${chrom[@]}; do
#       sample_gvcfs=""
#	for sample in $samples ; do
#	    sample_gvcfs=${sample_gvcfs}"-V $outdir/${sample}/gatk/${sample}.HC.${i}.g.vcf.gz \\"\n
#	done
#	time $gatk CombineGVCFs \
#        -R $reference/Homo_sapiens_assembly38.fasta \
#        ${sample_gvcgs} \
#        -O $outdir/population/${outname}.HC.${i}.g.vcf.gz && echo "** ${outname}.HC.${i}.g.vcf.gz done **" && \
#       time $gatk GenotypeGVCFs \
#        -R $reference/Homo_sapiens_assembly38.fasta \
#        -V $outdir/population/${outname.HC.${i}.g.vcf.gz \
#        -O $outdir/population/${outname}.HC.${i}.vcf.gz && echo "** ${outname}.HC.${i}.vcf.gz done **" &
#done && wait
#merge_vcfs=""
#for i in ${chrom[@]}; do
#       merge_vcfs=${merge_vcfs}" -I $outdir/population/${outname}.HC.${i}.vcf.gz \\"\n
#done && time $gatk MergeVcfs ${merge_vcfs} -O $outdir/population/${outname}.HC.vcf.gz && echo "** MergeVcfs done **"

#VQSR
#speed up for sample >500
#ignore sample information extract sites only
#time $gatk VariantRecalibrator \
#	-I $outdir/population/${outname}.HC.vcf.gz 、
#	-O $outdir/population/${outname}.HC.sites.vcf.gz && echo "** Only sites done **"

#SNP
time $gatk VariantRecalibrator \
        -R $reference/Homo_sapiens_assembly38.fasta \
        -V $outdir/population/${outname}.HC.vcf.gz \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf.gz \
        -resource:omini,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf.gz \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
        -mode SNP \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        --rscript-file $outdir/population/${outname}.HC.snps.plots.R \
        --tranches-fle $outdir/population/${outname}.HC.snps.tranches \
        -O $outdir/population/${outname}.HC.snps.recal && \
time $gatk ApplyVQSR \
        -R $reference/Homo_sapiens_assembly38.fasta \
        -V $outdir/population/${outname}.HC.vcf.gz \
        --ts_filter_level 99.0 \
        --tranches-file $outdir/population/${outname}.HC.snps.tranches \
        --recal-file $outdir/population/${outname}.HC.snps.recal \
	-mode SNP \
        -O $outdir/population/${outname}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"

#INDEL
time $gatk VariantRecalibrator \
        -R $reference/Homo_sapiens_assembly38.fasta \
        -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz \
        -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
        -mode INDEL \
        --max-gaussians 6 \
        --rscript-file $outdir/population/${outname}.HC.indels.plots.R \
        --tranches-fle $outdir/population/${outname}.HC.indels.tranches \
        -O $outdir/population/${outname}.HC.snps.indels.recal && \
time $gatk ApplyVQSR \
        -R $reference/Homo_sapiens_assembly38.fasta \
        -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz \
        --ts_filter_level 99.0 \
	--tranches-file $outdir/population/${outname}.HC.snps.indels.tranches \
        --recal-file $outdir/population/${outname}.HC.snps.indels.recal \
        -mode INDEL \
        -O $outdir/population/${outname}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR (${sample}.HC.VQSR.vcf.gz finish) done **"

#final results
#.g.vcf.gz;HC.vcf.gz;HC.VQSR.vcf.gz
