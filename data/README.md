Reference data

This directory should contain the ASCEND reference files required by the pipeline.

Large reference files are not stored directly in this GitHub repository. Download the ASCEND reference-data bundle from Zenodo:

TODO: add Zenodo DOI/URL

After downloading, unpack or copy the files into this directory:

tar -xzf ASCEND_reference_data_v1.tar.gz

The data/ directory should then contain:

dominant_genes_ENS.txt
gene_wFDR_weights_Sfacs50_quant.tsv
Uprod_dists.txt.gz
ALLVARS_MR_dist_by_pos5.txt.gz
ENS_ID2Gene_ID.txt.gz
BY_GENE_MR_5.txt.gz
muttargs5_noOL_RQC_format.txt.gz

These files provide the packaged mutation-rate annotations, gene-level mutation-rate expectations, gene-name mappings, missense mutation-rate coordinates, clustering reference distributions, and files used for weighted or censored false-discovery rate correction.

To check that the files are present, run from the main repository directory:

python src/check_setup.py
