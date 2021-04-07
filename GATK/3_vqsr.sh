#!/bin/bash
#fastq_to_gvcf

#path
gatk=/home/pang/program/gatk/gatk

#reference
#reference=/home/pang/data/public_data/GenomeRef/Human_ENSEMBLE_GRCh38/Sequence/GRCh38_r77.all.fa 
reference=/home/pang/program/gatk/GATK_bundle/Homo_sapiens_assembly38.fasta
GATK_bundle=/home/pang/program/gatk/GATK_bundle

## shell parameter
outname=$1 #prefix of out file

outdir=/home/pang/data/public_data/UKBB/exome_population

#population variation outdir
if [ ! -d $outdir/population ]
then mkdir -p $outdir/population
fi

#VQSR
#speed up for sample >500
#ignore sample information extract sites only
#time $gatk VariantRecalibrator \
#	-I $outdir/population/${outname}.HC.vcf.gz 、
#	-O $outdir/population/${outname}.HC.sites.vcf.gz && echo "** Only sites done **"

#SNP
time $gatk VariantRecalibrator \
        -R $reference \
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
        -R $reference \
        -V $outdir/population/${outname}.HC.vcf.gz \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file $outdir/population/${outname}.HC.snps.tranches \
        --recal-file $outdir/population/${outname}.HC.snps.recal \
	-mode SNP \
        -O $outdir/population/${outname}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"

#INDEL
time $gatk VariantRecalibrator \
        -R $reference \
        -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz \
        -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum \
        -mode INDEL \
        --max-gaussians 6 \
        --rscript-file $outdir/population/${outname}.HC.indels.plots.R \
        --tranches-fle $outdir/population/${outname}.HC.indels.tranches \
        -O $outdir/population/${outname}.HC.indels.recal && \
time $gatk ApplyVQSR \
        -R $reference \
        -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz \
        --truth-sensitivity-filter-level 99.0 \
	--tranches-file $outdir/population/${outname}.HC.indels.tranches \
        --recal-file $outdir/population/${outname}.HC.indels.recal \
        -mode INDEL \
        -O $outdir/population/${outname}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR (${sample}.HC.VQSR.vcf.gz finish) done **"

#final results
#.g.vcf.gz;HC.vcf.gz;HC.VQSR.vcf.gz
