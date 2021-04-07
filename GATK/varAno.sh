#!/bin/bash
#VEP annotation
VEP=/home/pang/program/ensembl-vep-release-98/vep
reference=/home/pang/program/gatk/GATK_bundle/Homo_sapiens_assembly38.fasta
#docker docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep
time $VEP --fasta $reference \
  	  --vcf --merged --fork 10 --hgvs --force_overwrite --everything \
	  --offline --dir_cache /home/pang/vep_data/homo_sapiens \
	  -i $outdir/gatk/${sample}.HC.VQSR.vcf.gz \
	  -o $outdir/gatk/${sample}.HC.VQSR.VEP.vcf.gz
