#!/bin/bash
#PDE5A located in chr4
#path
gatk=/home/pang/program/gatk/gatk

#reference
reference=/home/pang/data/public_data/GenomeRef/Human_ENSEMBLE_GRCh38/Sequence/GRCh38_r77.all.fa 
#reference=/home/pang/program/gatk/GATK_bundle/Homo_sapiens_assembly38.fasta
GATK_bundle=/home/pang/program/gatk/GATK_bundle

## shell parameter
sample_map=/home/pang/data/public_data/UKBB/ukbb_exome.sample_map
outname=$1 #prefix of out file
outdir=/home/pang/data/public_data/UKBB/exome_population

#population variation outdir
if [ ! -d $outdir/population ]
then mkdir -p $outdir/population
fi

time $gatk GenomicsDBImport \
        --sample-name-map ${sample_map} \
	--genomicsdb-workspace-path my_database \
	-L 4 \
	--batch-size 100 \
	--tmp-dir=/raid3/homedirs/pang/temp
time $gatk GenotypeGVCFs \
        -R $reference \
        -V gendb://$outdir/my_database \
        --use-new-qual-calculator \
        -O $outdir/population/${outname}.HC.vcf.gz && echo "** ${outname}.HC.vcf.gz done **"
