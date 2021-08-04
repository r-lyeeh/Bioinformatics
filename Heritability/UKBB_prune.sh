#UKBB is too large to calculate PCA or estimated grm, thus we used a pruned data set
1> Remove duplicate rs IDs
plink2 --bgen ukb_imp_chr*_QCed.bgen 'ref-first' --sample ukb_imp_chr0_QCed.sample --rm-dup 'force-first' --make-bed --out ukb_imp_chr*_QCed_rd
2> Pruned: optimizer refference from flashpca
plink2 --bfile ukb_imp_chr*_QCed_rd --indep-pairwise 1000 50 0.05 --exclude range /home/pang/program/flashpca/exclusion_regions_hg19.txt --out ukb_imp_chr*_QCed_rd
plink2 --bfile ukb_imp_chr*_QCed_rd --extract ukb_imp_chr*_QCed_rd.prune.in --make-bed --out ukb_imp_chr*
3> GRM extimated
gcta64 --bfile ukb_imp_chr* --make-grm-bin --out ukb_imp_chr*
4> Calculate PCA
gcta64 --mgrm multi_grm_ukb_imp.txt --make-grm-bin --out ukb_imp
gcta64 --grm ukb_imp --pca 20 --thread-num 4 --out ukb_imp_20PSs
