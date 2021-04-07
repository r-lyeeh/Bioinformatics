#!/bin/bash

#PBS -l nodes=ekgen6.isar-prestige.de:ppn=4
#PBS -q batch
#PBS -l meme=16gb

fimo_path=/raid3/homedirs/pang/program/meme-5.0.5/src/
excute_path=/raid3/homedirs/pang/data/MAFF/
cd ${excute_path}
${fimo_path}fimo -o MAFF_DNA_motif --verbosity 1 --thresh 1.0E-4 MAFF-LDLR-cofactor/MAFF_MA0495.1.meme MAFF-LDLR-cofactor/hg38.fa

# protein mast
${mast_path}mast -o MAFF_Protein_MAST -nostatus -minseqs 302 -remcorr -ev 10.0 $1 MAFF-LDLR-cofactor/Homo_sapiens.GRCh38.pep.all.fa

#protein meme
${meme_path}meme $1 -protein -o MAFF_Protein_MEME -nostatus -time 18000 -mod zoops -nmotifs 5 -minw 6 -maxw 50 -objfun classic -markov_order 0
