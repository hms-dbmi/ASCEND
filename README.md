# ASCEND

ASCEND is a Python tool for gene-level association testing from observed *de novo* variation. The pipeline converts annotated variant-level observations into cohort-level summary statistics, adds packaged mutation-rate reference annotations, computes component association statistics, combines evidence across variant classes, and performs weighted false-discovery rate control.

## Repository structure

```text
ASCEND_tool/
  README.md
  LICENSE
  requirements.txt
  .gitignore

  data/
    dominant_genes_ENS.txt
    gene_wFDR_weights_Sfacs50_quant.tsv
    Uprod_dists.txt.gz
    ALLVARS_MR_dist_by_pos5.txt.gz
    ENS_ID2Gene_ID.txt.gz
    BY_GENE_MR_5.txt.gz
    muttargs5_noOL_RQC_format.txt.gz

  src/
    ASCEND_main.py
    preprocess.py
    reference.py
    statistics.py
    combine.py
    fdr.py
    check_setup.py

  example/
    test_vars.vcf.gz
    lazy-test.sh
```

## Installation

ASCEND requires Python 3.10 or later and the following Python packages:

- NumPy
- SciPy
- pandas

A minimal conda environment can be created with:

```bash
conda create -n ascend python=3.10 numpy scipy pandas
conda activate ascend
```

Alternatively, install the Python dependencies with:

```bash
pip install -r requirements.txt
```

The analyses described in the manuscript were run using Python 3.10.18 with NumPy 2.2.6, SciPy 1.14.0, and pandas 2.3.2.

## Reference data

ASCEND requires packaged reference files in the `data/` directory. These files include variant annotations, gene-level mutation-rate expectations, per-position mutation-rate coordinates, gene-symbol mappings, and files used for weighted or censored false-discovery rate control.

Expected files:

```text
dominant_genes_ENS.txt
gene_wFDR_weights_Sfacs50_quant.tsv
Uprod_dists.txt.gz
ALLVARS_MR_dist_by_pos5.txt.gz
ENS_ID2Gene_ID.txt.gz
BY_GENE_MR_5.txt.gz
muttargs5_noOL_RQC_format.txt.gz
```

Large reference files are not intended to be stored directly in the GitHub repository. Download the ASCEND reference-data bundle from:

```text
TODO: add Zenodo/Figshare/Dataverse DOI or URL
```

Then unpack or copy the files into:

```text
ASCEND_tool/data/
```

After adding the reference files, check the installation with:

```bash
python ASCEND_tool/check_setup.py
```

## Basic usage

Run ASCEND from a VCF-like file containing observed *de novo* variants:

```bash
python src/ASCEND_main.py --i path/to/cohort.vcf.gz --o cohort_out
```

Run ASCEND from precomputed per-gene summary statistics:

```bash
python src/ASCEND_main.py --i cohort_out_sumstats.txt --o cohort_from_sumstats
```

By default, ASCEND expects reference files in `../data` relative to `src/ASCEND_main.py`. If needed, use the command-line options in `ASCEND_main.py --help` to specify alternative reference-file locations.

## Input files

ASCEND can be run either from a VCF-like file or from precomputed per-gene summary statistics.

For VCF-like input, each non-header row is treated as one observed *de novo* variant. The file must contain chromosome, position, reference allele, and alternate allele fields in standard VCF-like columns. Variants are matched to the packaged ASCEND annotation table using the chromosome-position-reference-alternate key. Recurrent observations of the same variant are counted as separate observations.

For summary-statistics input, the file should contain one row per gene with observed loss-of-function and synonymous variant counts, missense score sums, and observed missense variant positions. The expected core columns are:

```text
ENS_ID
AM_y
REVEL_y
PAI_y
Lof_varN
syn_varN
missense_positions
```

Missing values should be encoded as `NA`.

## Output files

A standard ASCEND run writes files using the prefix supplied with `--o`. Typical outputs include:

```text
<prefix>_sumstats.txt
<prefix>_sumstats_results.txt
```

Depending on input type and command-line options, intermediate files may also be written, including observed variant tables, reference-annotated summary statistics, and component P-value files.

## Example

A small toy example is provided in the `example/` directory. test by running:

```bash
cd example
bash lazy-test.sh
```

This script tests ASCEND from VCF input, from precomputed summary statistics, and with alternative clustering approximations.
