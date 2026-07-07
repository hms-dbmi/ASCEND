# ASCEND reference data

This directory is intentionally incomplete in the GitHub repository.

The large ASCEND reference files required by the pipeline are **not stored on GitHub**. They must be downloaded from Zenodo:

```text
https://doi.org/10.5281/zenodo.21249781
```

The Zenodo record contains the file:

```text
ASCEND_data_dir.zip
```

Download this archive and unpack it into the main ASCEND repository directory:

```bash
cd ASCEND
unzip ASCEND_data_dir.zip
```

After unpacking, this directory should contain:

```text
dominant_genes_ENS.txt
gene_wFDR_weights_Sfacs50_quant.tsv
Uprod_dists.txt.gz
ALLVARS_MR_dist_by_pos5.txt.gz
ENS_ID2Gene_ID.txt.gz
BY_GENE_MR_5.txt.gz
muttargs5_noOL_RQC_format.txt.gz
```

These files provide the packaged mutation-rate annotations, gene-level mutation-rate expectations, gene-name mappings, mutation-rate-scaled missense coordinates, clustering reference distributions, and files used for weighted or censored false-discovery rate correction.

To check that all required files are present, run from the repository root:

```bash
python check_setup.py
```

The large reference files should not be committed to GitHub.
