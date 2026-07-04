# ASCEND source code

This directory contains the main ASCEND source code.

## Main entry point

```bash
python ASCEND_main.py --i <input_file> --o <output_prefix>
```

`ASCEND_main.py` is the main command-line interface. It can run ASCEND either from a VCF-like file containing observed *de novo* variants or from a precomputed per-gene summary-statistics file.

By default, the script expects ASCEND reference files to be located in `../data/` relative to this directory.

## Source files

* `ASCEND_main.py`
  Main pipeline script. Coordinates preprocessing, reference annotation, association testing, P-value combination, and false-discovery rate correction.

* `preprocess.py`
  Converts VCF-like variant input into observed per-gene summary statistics.

* `reference.py`
  Adds packaged ASCEND reference annotations, including mutation-rate expectations, gene names, and mutation-rate-scaled missense coordinates.

* `statistics.py`
  Computes component gene-level association P-values, including loss-of-function enrichment, missense enrichment, and missense clustering.

* `combine.py`
  Combines component P-values across missense scores and variant classes.

* `fdr.py`
  Applies weighted and censored false-discovery rate correction.

* `check_setup.py`
  Checks that required Python packages and ASCEND reference files are available.

## Example

From the `example/` directory:

```bash
bash lazy-test.sh
```

or directly from this directory:

```bash
python ASCEND_main.py --i ../example/test_vars.vcf.gz --o ../example/test_vars_out
```
