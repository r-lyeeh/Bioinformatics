#!/bin/bash
set=$1
perl convert2annovar.pl -format vcf4 -allsample ${set}.recode.vcf -withfreq >${set}.avinput_raw
perl table_annovar.pl ${set}.avinput ~/program/annovar/humandb/ -buildver hg19 -out ${set} -remove -protocol refGene -operation g -nastring . -csvout -polish
