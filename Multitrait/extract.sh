#!/bin/bash
cd /home/pang/Desktop/CAD_IS_PAD_20200121/data/SNP/5000
for i in cad is pad
do
awk '{print $3}' ${i}_kb500_gene >${i}_raw_snp
awk '{print $16}' ${i}_kb500_gene >${i}_gene_raw
sort ${i}_gene_raw >t1
uniq t1 >${i}_gene_raw
awk '{print $14}' ${i}_locus >${i}_locus_cytoband
sed -i -e 's/cad//g' ${i}_locus_cytoband
sed -i -e 's/is//g' ${i}_locus_cytoband
sed -i -e 's/pad//g' ${i}_locus_cytoband
sort ${i}_locus_cytoband >t1
uniq t1 >${i}_locus_cytoband
awk '{print $14}' ${i}_locus_1m >${i}_locus_1m_raw
sed -i -e 's/cad//g' ${i}_locus_1m_raw
sed -i -e 's/is//g' ${i}_locus_1m_raw
sed -i -e 's/pad//g' ${i}_locus_1m_raw
sort ${i}_locus_1m_raw >t1
uniq t1 >${i}_locus_1m_raw
done
grep -Fxf cad_raw_snp is_raw_snp >cad_is_snp
grep -Fxf cad_raw_snp pad_raw_snp >cad_pad_snp
grep -Fxf is_raw_snp pad_raw_snp >is_pad_snp
grep -Fxf cad_is_snp pad_raw_snp >cad_is_pad_snp
grep -Fxf cad_locus_cytoband is_locus_cytoband >cad_is_locus_cytoband
grep -Fxf cad_locus_cytoband pad_locus_cytoband >cad_pad_locus_cytoband
grep -Fxf is_locus_cytoband pad_locus_cytoband >is_pad_locus_cytoband
grep -Fxf cad_is_locus_cytoband pad_locus_cytoband >cad_is_pad_locus_cytoband
grep -Fxf cad_locus_1m_raw is_locus_1m_raw >cad_is_locus_1m_raw
grep -Fxf cad_locus_1m_raw pad_locus_1m_raw >cad_pad_locus_1m_raw
grep -Fxf is_locus_1m_raw pad_locus_1m_raw >is_pad_locus_1m_raw
grep -Fxf cad_is_locus_1m_raw pad_locus_1m_raw >cad_is_pad_locus_1m_raw
grep -Fxf cad_gene_raw is_gene_raw >cad_is_gene_raw
grep -Fxf cad_gene_raw pad_gene_raw >cad_pad_gene_raw
grep -Fxf is_gene_raw pad_gene_raw >is_pad_gene_raw
grep -Fxf cad_is_gene_raw pad_gene_raw >cad_is_pad_gene_raw
