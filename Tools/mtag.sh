source activate mtag
python2 ../mtag.py \
        --sumstats 1_OA2016_hm3samp_NEUR.txt,1_OA2016_hm3samp_SWB.txt \
        --out ./tutorial_results_1.1NS \
        --n_min 0.0 \
        --use_beta_se \
        --stream_stdout &
        # --no_overlap # When there is no sample overlap in the GWAS
        # --equal_h2 # Performing meta-analysis with mtag # assume variation between traits is only due to non-genetic factors
# maxFDR calculation
# --fdr \
# --stream_stdout

### Necessary columns for summary statistics
# rs number; number of studies; p-value; Allele 1; Allele 2; N; Sgned statistics (z, or, beata, se); info; maf
