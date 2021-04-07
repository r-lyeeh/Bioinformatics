#~/bin/bash
#PBS -q batch
#PBS -j oe
cd /home/pang/data/public_data/UKBB/exome
list=`ls *_23161_0_0.gvcf.gz`
for file in ${list}
do
#	sample=${file%_23161_0_0.gvcf.gz}
	tabix -p vcf $file
#	zcat $file |bgzip -c >new.$sample.vcf.gz && tabix new.$sample.vcf.gz #if above failed
done
