# ASCEND-Diag

ASCEND-Diag estimates diagnostic probabilities for observed *de novo* variants using a reference disease cohort, mutation-rate annotations, variant functional class, and missense pathogenicity scores.

ASCEND-Diag is distinct from the main ASCEND association-testing workflow. The main ASCEND pipeline tests whether genes are enriched for *de novo* variation. ASCEND-Diag instead estimates the probability that a specific observed variant contributes to the phenotype of an individual proband, under the assumptions of the ASCEND-Diag model.

## Script

The main script is:

```bash
ASCEND_diag.py
```

## Required inputs

ASCEND-Diag requires four main inputs:

### 1. Gene list

A text file containing the genes to evaluate, with one Ensembl gene ID per line.

Example:

```text
ENSG00000141510
ENSG00000155657
ENSG00000171862
```

Lines beginning with `#` are ignored.

### 2. Reference cohort variant-info file

A tab-separated file containing annotated variants observed in the reference disease cohort. This is typically generated during ASCEND preprocessing.

Required columns include:

```text
ensembl_gene_id
SAMPLE_RECURRENCE
Synonymous-flag
LOFTEE-flag
AlphaMissense_MT
REVEL_MT
PrimateAI-3D_MT
```

### 3. Reference cohort per-gene summary-statistics file

A tab-separated file containing per-gene summary statistics and mutation-rate expectations for the reference disease cohort.

Required columns include:

```text
ENS_ID
syn_varN
Lof_varN
Lof_lambda
AM_lambda
REVEL_lambda
PAI_lambda
```

### 4. Processed target-variant file

A tab-separated file containing the variants for which diagnostic probabilities should be computed. This file may be gzipped.

Required columns include:

```text
ensembl_gene_id
Synonymous-flag
LOFTEE-flag
AlphaMissense_MT
REVEL_MT
PrimateAI-3D_MT
```

Additional columns are preserved in the output.

## Basic usage

```bash
python ASCEND_diag.py \
  --gene-list-file diagnostic_genes.txt \
  --ref-var-info-file reference_variant_info.tsv \
  --ref-sum-stats-file reference_sumstats.tsv \
  --processed-vars-file target_variants.tsv.gz \
  --output-file target_variants_diag.tsv.gz \
  --colname DiagProb \
  --P_D 0.01 \
  --pnull 0.25
```

## Command-line options

```text
--gene-list-file
    File containing the gene whitelist, one Ensembl gene ID per line.
    Default: gene_list_file.txt

--ref-var-info-file
    Reference cohort variant-info file.

--ref-sum-stats-file
    Reference cohort per-gene summary-statistics file.

--processed-vars-file
    Variants for which diagnostic probabilities should be calculated.
    May be plain text or gzip-compressed.

--output-file
    Output file. If the filename ends in .gz, output is gzip-compressed.

--colname
    Prefix for added diagnostic-probability columns.
    Default: DiagProb

--P_D
    Assumed prevalence of the umbrella phenotype.
    Default: 0.01

--pnull
    Assumed probability of LoF misannotation.
    Default: 0.25

--add-one-to-varN
    Add one to the observed variant count in enrichment calculations.

--verbose
    Print progress messages.
```

## Output

ASCEND-Diag writes the input target-variant table with additional diagnostic-probability columns appended.

For each missense score, four diagnostic-probability estimates are reported:

```text
<colname>_PP_AlphaMissense_MT_1
<colname>_PP_AlphaMissense_MT_2
<colname>_PP_AlphaMissense_MT_3
<colname>_PP_AlphaMissense_MT_4

<colname>_PP_REVEL_MT_1
<colname>_PP_REVEL_MT_2
<colname>_PP_REVEL_MT_3
<colname>_PP_REVEL_MT_4

<colname>_PP_PrimateAI-3D_MT_1
<colname>_PP_PrimateAI-3D_MT_2
<colname>_PP_PrimateAI-3D_MT_3
<colname>_PP_PrimateAI-3D_MT_4
```

For LoF variants, the enrichment-based LoF diagnostic probability is copied across these output columns. For synonymous variants and variants that are neither LoF nor scored missense, diagnostic-probability columns are set to `NA`.

The four missense estimates correspond to different ways of estimating diagnostic probability:

```text
1. Gene-level missense enrichment.
2. Genome-wide score-based pathogenicity model.
3. Gene-specific score-based pathogenicity model.
4. Combined score and enrichment-based model.
```

## Reference data

ASCEND-Diag uses the ASCEND mutation-rate and annotation framework. Some required annotation files are shared with the main ASCEND pipeline, including:

```text
muttargs5_noOL_RQC_format.txt.gz
```

Additional diagnostic-reference files may be required, depending on the reference cohort used. These files are not necessarily included in the public GitHub repository, especially if they are derived from restricted-access clinical or cohort data.

If public diagnostic reference files are available, download them from:

```text
TODO: add Zenodo DOI/URL
```

and place them in the appropriate data directory.

## Notes and limitations

ASCEND-Diag assumes that the incoming proband and the reference disease cohort are ascertained under a comparable umbrella phenotype. Diagnostic probabilities therefore depend on the choice of reference cohort, the assumed phenotype prevalence `P_D`, the variant annotation model, and the missense pathogenicity score calibration.

ASCEND-Diag is research-use software. Its output should not be used as the sole basis for clinical diagnosis or medical decision-making.
