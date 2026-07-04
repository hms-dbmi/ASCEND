This directory contains a small toy example for testing the main ASCEND workflow.

Files
test_vars.vcf.gz
Toy VCF-like file containing observed de novo variants.
lazy-test.sh
Test script that runs ASCEND in several modes:
from VCF-like variant input;
from precomputed per-gene summary statistics;
using the chi-square clustering approximation;
using the permutation-based clustering approximation.
Running the example

From this directory, run:

bash lazy-test.sh

The script calls the main ASCEND pipeline in ../src/ASCEND_main.py
