#!/usr/bin/env bash
set -euo pipefail

# Run from the directory containing this script
cd "$(dirname "$0")"

echo "Testing ASCEND from VCF input..."
python ../src/ASCEND_main.py --i=test_vars.vcf.gz --o=test_vars_out

echo "Testing ASCEND from precomputed summary statistics..."
python ../src/ASCEND_main.py --i=test_vars_out_sumstats.txt --o=test_ss_out

echo "Testing chi-square clustering approximation..."
python ../src/ASCEND_main.py --i=test_vars_out_sumstats.txt --o=test_chi2_out --clust_method=chi2

echo "Testing permutation clustering approximation..."
python ../src/ASCEND_main.py --i=test_vars_out_sumstats.txt --o=test_permut_out --clust_method=permut

echo "Done. Generated files:"
ls -lh test_vars_out* test_ss_out* test_chi2_out* test_permut_out*
