#!/bin/bash
#PDE5A located in chr4
#path
gatk=/home/pang/program/gatk/gatk

#reference
reference=/home/pang/program/gatk/GATK_bundle/Homo_sapiens_assembly38.fasta
GATK_bundle=/home/pang/program/gatk/GATK_bundle

## shell parameter
sample_map=/home/pang/data/public_data/UKBB/ukbb_exome.sample_map
outname=$1 #prefix of out file
outdir=/home/pang/data/SCP_involve/PDE5A_Thorsten_20191011/gatk

#population variation outdir
if [ ! -d $outdir/population ]
then mkdir -p $outdir/population
fi

#Merge samples using sample_map
time $gatk GenomicsDBImport \
        --sample-name-map ${sample_map} \
        --genomicsdb-workspace-path my_database \
	--batch-size 100
	--intervals chr:ds1-ds2 \
	--tmp-dir=/raid3/homedirs/pang/temp
time $gatk GenotypeGVCFs \
	-R $reference \
	-V $outdir/population/mydatabase \
	--use-new-qual-calculator \
	-O $outdir/population/${outname}.HC.vcf.gz && echo "** ${outname}.HC.vcf.gz done **"
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
        --tranches-file $outdir/population/${outname}.HC.snps.tranches \
        -O $outdir/population/${outname}.HC.snps.recal && \
time $gatk ApplyVQSR \
        -R $reference \
        -V $outdir/population/${outname}.HC.vcf.gz \
        --ts_filter_level 99.0 \
        --tranches-file $outdir/population/${outname}.HC.snps.tranches \
        --recal-file $outdir/population/${outname}.HC.snps.recal \
        -mode SNP \
        -O $outdir/population/${outname}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"

#INDEL
time $gatk VariantRecalibrator \
        -R $reference \
        -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz \
        -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
        -mode INDEL \
        --max-gaussians 6 \
        --rscript-file $outdir/population/${outname}.HC.indels.plots.R \
        --tranches-file $outdir/population/${outname}.HC.indels.tranches \
        -O $outdir/population/${outname}.HC.snps.indels.recal && \
time $gatk ApplyVQSR \
        -R $reference \
        -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz \
        --ts_filter_level 99.0 \
        --tranches-file $outdir/population/${outname}.HC.snps.indels.tranches \
        --recal-file $outdir/population/${outname}.HC.snps.indels.recal \
        -mode INDEL \
        -O $outdir/population/${outname}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR (${sample}.HC.VQSR.vcf.gz finish) done **"

