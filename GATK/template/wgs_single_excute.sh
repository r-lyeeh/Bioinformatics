#!/bin/bash
## single sampe Illunina, PE fasta

#path
trimmomatic=/home/pang/program/Trimmomatic-0.39/trimmomatic-0.39.jar
bwa=/home/pang/program/bwa/bwa
samtools=/home/pang/program/samtools/samtools
gatk=/home/pang/program/gatk/gatk

#reference
reference=/home/pang/program/gatk/GATK_bundle/Homo_sapiens_assembly38.fasta
GATK_bundle=/home/pang/program/gatk/GATK_bundle

## shell parameter
fq1=$1
fq2=$2
RGID=$3 #Read Group ~ Lane ID
library=$4 #Library number
sample=$5 #sample ID
outdir=$6 #output path

##Sample outdir
outdir=${outdir}/${sample}

##fastq name
fq_file_name=`basename $fq1`
fq_file_name=${fq1_file_name%%.1.fq.gz}

# output directory
if [ ! -d $outdir/cleanfq ]
then mkdir -p $outdir/cleanfq
fi

if [ ! -d $outdir/bwa ]
then mkdir -p $outdir/bwa
fi

if [ ! -d $outdir/gatk ]
then mkdir -p $outdir/gatk
fi

##Trimmomatic QC keepBothReads=True
time java -jar ${trimmomatic} PE \
	$fq1 $fq2 \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz ${fq_file_name}.unpaired.1.fq.gz \
	$outdir/cleanfq/${fq_file_name}.paired.2.fq.gz ${fq_file_name}.unpaired.2.fq.gz \
	ILLUMINACLIP:/home/pang/program/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:8:True \
	SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 && echo "** fq QC done **"

##bwa mem alignment
time $bwa mem -t 8 -M -Y -R "@RG\tID:$RGID\tPL:ILLUMINA\tPU:$PU\tLB:$library\tSM:$sample" $reference/Homo_sapiens_assembly38.fasta \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz $outdir/cleanfq/${fq_file_name}.paired.2.fq.gz | $samtools view -Sb - > $outdir/bwa/${sample}.bam && \
	echo "** BWA MEM done **" && \
time $samtools sort -@ 4 -m 4G -O bam -o $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.bam && echo "** sorted raw bamfile done **"
# time $samtools index $ourdir/bwa/${sample}.sorted.bam && echo "** ${sample}.sorted.bam index done **"

## mark duplicate seq
$gatk MarkDuplicates \
	-I $outdir/bwa/${sample}.sorted.bam \
	-M $outdir/bwa/${sample}.markdup_metrics.txt \
	-O $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.bam MarkDuplicates done **"
## Index for markdup
time $samtools index $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.markdup.bam index done **"

## BQSR INDEX necessary
time $gatk BaseRecalibrator \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	--known-sites $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--known-sites $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
	-O $outdir/bwa/${sample}.sorted.markdup.recal_data.table && echo "** ${sample}.sorted.markdup.recal_data.table done **"

time $gatk ApplyBQSR \
	--bqsr-recal-file $outdir/bwa/${sample}.sorted.markdup.recal_data.table \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	-O $outdir/bwa/${sample}.sorted.markdip.BQSR.bam && echo "** ApplyBQSR done **"

time $samtools index $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && echo "** ${sample}.sorted.markdup.BQSR.bam index done **"

## detect variation
time $gatk HaplotypeCaller \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.vcf.gz && echo "** ${sample}.HC.vcf.gz done **"

## sampel chr -> joint chr
#chrom=(chr1 chr2 chr3 chr4...
#for i in ${chrom[@]}; do
#	time $gatk HaplotypeCaller \
#        -R $reference/Homo_sapiens_assembly38.fasta \
#        -I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#	 -L $i \
#        -O $outdir/gatk/${sample}.HC.${i}.vcf.gz && echo "** ${sample}.HC.${i}.vcf.gz done ** &
#done && wait
#merge_vcfs=""
#for i in ${chrom[@]}; do
#	merge_vcfs=${merge_vcfs}" -I $outdir/gatk/${sample}.HC.${i}.vcf.gz \\"\n
#done && time $gatk MergeVcfs ${merge_vcfs} -O $outdir/gatk/${sample}.HC.vcf.gz && echo "** MergeVcfs done **"

#Population gVCF -> GenotypeGVCFs multi-sample VCF
#time $gatk HaplotypeCaller \
#	--emit-ref-confidence GVCF \
#        -R $reference/Homo_sapiens_assembly38.fasta \
#        -I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#        -O $outdir/gatk/${sample}.HC.g.vcf.gz && echo "** GVCF ${sample}.HC.g.vcf.gz done **" 
#time $gatk GenotypeGVCFs \
#        -R $reference/Homo_sapiens_assembly38.fasta \
#	 -V $outdir/gatk/${sample}.HC.g.vcf.gz \
#	 -O $outdir/gatk/${sample}.HC.vcf.gz && echo "** ${sample}.HC.vcf.gz done **"

#chrom vcf - genotype gvcf
#chrom=(chr1 chr2 chr3 chr4...
#for i in ${chrom[@]}; do
#       time $gatk HaplotypeCaller \
#        -R $reference/Homo_sapiens_assembly38.fasta \
#        -I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#        -L $i \
#        -O $outdir/gatk/${sample}.HC.${i}.vcf.gz && \
#	time $gatk GenotypeGVCFs \
#        -R $reference/Homo_sapiens_assembly38.fasta \
#        -V $outdir/gatk/${sample}.HC.${i}.g.vcf.gz \
#        -O $outdir/gatk/${sample}.HC.${i}.vcf.gz && echo "** ${sample}.HC.${i}.vcf.gz done **" &
#done && wait
#merge_vcfs=""
#for i in ${chrom[@]}; do
#       merge_vcfs=${merge_vcfs}" -I $outdir/gatk/${sample}.HC.${i}.vcf.gz \\"\n
#done && time $gatk MergeVcfs ${merge_vcfs} -O $outdir/gatk/${sample}.HC.vcf.gz && echo "** MergeVcfs done **"

#VQSR SNP&INDEL has different standard
#SNP
time $gatk VariantRecalibrator \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf.gz \
	-resource:omini,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf.gz \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--rscript-file $outdir/gatk/${sample}.HC.snps.plots.R \
	--tranches-fle $outdir/gatk/${sample}.HC.snps.tranches \
	-O $outdir/gatk/${sample}.HC.snps.recal && \
time $gatk ApplyVQSR \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	--ts_filter_level 99.0 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	--recal-file $outdir/gatk/${sample}.HC.snps.recal \
	-mode SNP \
	-O $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"

#INDEL
time $gatk VariantRecalibrator \
        -R $reference/Homo_sapiens_assembly38.fasta \
        -V $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
        -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
        -mode INDEL \
	--max-gaussians 6 \
        --rscript-file $outdir/gatk/${sample}.HC.indels.plots.R \
        --tranches-fle $outdir/gatk/${sample}.HC.indels.tranches \
        -O $outdir/gatk/${sample}.HC.snps.indels.recal && \
time $gatk ApplyVQSR \
        -R $reference/Homo_sapiens_assembly38.fasta \
        -V $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
        --ts_filter_level 99.0 \
        --tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
        --recal-file $outdir/gatk/${sample}.HC.snps.indels.recal \
        -mode INDEL \
        -O $outdir/gatk/${sample}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR (${sample}.HC.VQSR.vcf.gz finish) done **"

#final results
#${sample}.sorted.markdup.BQSR.bam
#.g.vcf.gz; HC.vcf.gz; HC.VQSR.vcg.gz

