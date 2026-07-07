# ASCEND

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21249781.svg)](https://doi.org/10.5281/zenodo.21249781)

ASCEND (Association Statistics for Clustering and ENrichment of De novo variants) is a Python tool for gene-level association testing from observed *de novo* variation. The pipeline converts annotated variant-level observations into cohort-level summary statistics, adds packaged mutation-rate reference annotations, computes component association statistics, combines evidence across variant classes, and performs weighted false-discovery rate control.

## Repository structure

```text
ASCEND/
  README.md
  requirements.txt
  check_setup.py

  data/
    README.md
    # large reference files downloaded from Zenodo

  src/
    ASCEND_main.py
    preprocess.py
    reference.py
    statistics.py
    combine.py
    fdr.py

  diag/
    ASCEND_Diag.py
    README.md
    # ASCEND-Diag reference cohort files

  example/
    test_vars.vcf.gz
    lazy-test.sh
    README.md
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

ASCEND requires large packaged reference files that are **not stored in this GitHub repository**. The GitHub repository contains the source code, documentation, and toy example only. To run ASCEND on real data, download the ASCEND data archive from Zenodo:

```text
https://doi.org/10.5281/zenodo.21249781
```

The Zenodo record contains `ASCEND_data_dir.zip`, which includes the complete `data/` directory required by the ASCEND pipeline.

After cloning the repository, download and unpack the data archive into the repository root:

```bash
git clone https://github.com/hms-dbmi/ASCEND.git
cd ASCEND

# Download ASCEND_data_dir.zip from Zenodo, then run:
unzip ASCEND_data_dir.zip
```

After unpacking, the repository should contain:

```text
ASCEND/
  data/
    dominant_genes_ENS.txt
    gene_wFDR_weights_Sfacs50_quant.tsv
    Uprod_dists.txt.gz
    ALLVARS_MR_dist_by_pos5.txt.gz
    ENS_ID2Gene_ID.txt.gz
    BY_GENE_MR_5.txt.gz
    muttargs5_noOL_RQC_format.txt.gz
```

Check that the reference files and Python dependencies are available with:

```bash
python check_setup.py
```

If the setup check succeeds, the toy example can be run with:

```bash
cd example
bash lazy-test.sh
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

## Data availability

The ASCEND Zenodo record contains:

- `ASCEND_data_dir.zip`: reference files required to run ASCEND;
- `GeneInfo.tsv.gz`: precomputed gene-level annotations and ASCEND association statistics;
- `var_diag_probs.txt.gz`: precomputed ASCEND-Diag diagnostic probabilities for coding variants;
- file-specific README files describing these tables.

Zenodo DOI: https://doi.org/10.5281/zenodo.21249781

This script tests ASCEND from VCF input, from precomputed summary statistics, and with alternative clustering approximations.
