"""
pipeline.py

End-to-end ASCEND workflow.

This module wires together:

    preprocess.py   : VCF -> observed per-gene summary stats
    reference.py    : add mutation-rate reference, gene names, MR coordinates
    statistics.py   : add component association P-values
    combine.py      : combine component P-values, if available
    fdr.py          : weighted/censored FDR, if requested

    ../data/

relative to the location of this pipeline.py file.

All default paths can still be overridden from the command line.
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Default data/config paths
# ---------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = (SCRIPT_DIR / "../data").resolve()

DEFAULT_DATA_FILES = {
    "annotation_file": "muttargs5_noOL_RQC_format.txt.gz",
    "by_gene_MR_file": "BY_GENE_MR_5.txt.gz",
    "ENS2GeneID_file": "ENS_ID2Gene_ID.txt.gz",
    "MR_dist_file": "ALLVARS_MR_dist_by_pos5.txt.gz",
    "Uprod_file": "Uprod_dists.txt.gz",
    "gene_constraint_file": "gene_wFDR_weights_Sfacs50_quant.tsv",
    "OMIM_genes_file": "dominant_genes_ENS.txt",
}


def resolve_data_paths(data_dir: str | None = None) -> dict[str, str]:
    """
    Return ASCEND packaged-data paths.

    If data_dir is omitted, paths are resolved as ../data relative to this file.
    """
    base = Path(data_dir).resolve() if data_dir is not None else DEFAULT_DATA_DIR
    return {key: str(base / filename) for key, filename in DEFAULT_DATA_FILES.items()}


def use_default_if_none(value, default_value):
    """Return default_value only when value is None."""
    return default_value if value is None else value


# ---------------------------------------------------------------------
# Imports of ASCEND modules
# ---------------------------------------------------------------------

try:
    from . import preprocess
    from . import reference
    from . import statistics
except ImportError:
    import preprocess
    import reference
    import statistics


try:
    from . import fdr
except ImportError:
    try:
        import fdr
    except ImportError:
        fdr = None


def write_tsv(df: pd.DataFrame, output_file: str) -> None:
    """Write a DataFrame as TSV using ASCEND's NA convention."""
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, sep="\t", index=False, na_rep="NA")


def read_tsv(path_str: str) -> pd.DataFrame:
    """Read TSV/TSV.GZ with compression inferred from filename."""
    return pd.read_table(path_str, sep="\t", compression="infer")


def validate_columns(columns, required, source_name: str = "") -> None:
    """Raise a readable error if required columns are missing."""
    missing = set(required) - set(columns)
    if missing:
        label = f" in {source_name}" if source_name else ""
        raise ValueError(f"Missing required columns{label}: {sorted(missing)}")


def require_existing_file(path_str: str | None, label: str) -> None:
    """
    Fail early with a readable message if a required reference/config file is absent.
    """
    if path_str is None:
        raise ValueError(f"{label} is required but was not supplied.")

    if not Path(path_str).exists():
        raise FileNotFoundError(
            f"{label} does not exist: {path_str}\n"
            "Either place the file in ../data relative to pipeline.py, pass --data_dir, "
            "or override this path explicitly."
        )


def split_input_files(input_files):
    """
    Split input files into VCFs and non-VCF summary-stat files.
    """
    sumstat_files = []
    vcf_files = []

    for file_name in input_files:
        file_name = file_name.strip()
        if not file_name:
            continue

        if file_name.endswith(".vcf") or file_name.endswith(".vcf.gz"):
            vcf_files.append(file_name)
        else:
            sumstat_files.append(file_name)

    return {"vcfs": vcf_files, "sumstats": sumstat_files}


def file_mask_from_path(path_str: str) -> str:
    """Create a stable output mask from an input filename."""
    base = os.path.basename(path_str)
    for suffix in [".vcf.gz", ".vcf", ".tsv.gz", ".tsv", ".txt.gz", ".txt"]:
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return base


def require_reference_args(
    by_gene_MR_file: str | None,
    ENS2GeneID_file: str | None,
    MR_dist_file: str | None,
) -> None:
    """Check that reference files were supplied and exist."""
    require_existing_file(by_gene_MR_file, "by-gene mutation-rate reference file")
    require_existing_file(ENS2GeneID_file, "ENS_ID-to-gene-name mapping file")
    require_existing_file(MR_dist_file, "per-position MR-distance reference file")


def preprocess_one_vcf(
    vcf_file: str,
    annotation_file: str,
    out_prefix: str,
    write_variant_table: bool = True,
) -> str:
    """
    Run preprocess.py on one VCF and return the observed summary-stat filename.
    """
    require_existing_file(annotation_file, "variant annotation/mutation-target file")

    vcf_mask = file_mask_from_path(vcf_file)
    sample_prefix = f"{out_prefix}_{vcf_mask}"

    gene_stats, observed_variants = preprocess.preprocess_vcf(
        vcf_file=vcf_file,
        annotation_file=annotation_file,
        keep_variant_table=write_variant_table,
    )

    observed_sumstats_file = f"{sample_prefix}_observed_sumstats.tsv"
    write_tsv(gene_stats, observed_sumstats_file)

    if write_variant_table:
        observed_variants_file = f"{sample_prefix}_observed_variants.tsv"
        write_tsv(observed_variants, observed_variants_file)

    return observed_sumstats_file


def add_reference_to_observed_sumstats(
    observed_sumstats_file: str,
    by_gene_MR_file: str,
    ENS2GeneID_file: str,
    MR_dist_file: str,
    out_prefix: str,
) -> str:
    """
    Run reference.py on one observed summary-stat file and return the output filename.
    """
    require_reference_args(by_gene_MR_file, ENS2GeneID_file, MR_dist_file)

    sample_mask = file_mask_from_path(observed_sumstats_file)
    reference_file = f"{out_prefix}_{sample_mask}_reference_annotated.tsv"

    if hasattr(reference, "update_sumstats_file_with_reference"):
        reference.update_sumstats_file_with_reference(
            observed_sumstats_file=observed_sumstats_file,
            by_gene_mr_file=by_gene_MR_file,
            gene_name_file=ENS2GeneID_file,
            mr_dist_file=MR_dist_file,
            output_file=reference_file,
        )
    elif hasattr(reference, "update_sum_stats"):
        reference.update_sum_stats(
            observed_sumstats_file,
            by_gene_MR_file,
            ENS2GeneID_file,
            MR_dist_file,
            reference_file,
        )
    else:
        raise AttributeError(
            "reference.py must define update_sumstats_file_with_reference(...) "
            "or update_sum_stats(...)."
        )

    return reference_file


def detect_sumstats_stage(sumstats_file: str) -> str:
    """
    Guess whether a summary-stat file is already reference-annotated.

    Returns "reference" or "observed".
    """
    df_head = pd.read_table(sumstats_file, sep="\t", compression="infer", nrows=5)
    cols = set(df_head.columns)

    reference_markers = {
        "AM_lambda",
        "REVEL_lambda",
        "PAI_lambda",
        "Lof_lambda",
        "syn_lambda",
        "missense_lambda",
        "missense_MTs",
        "missense_MTs_by_pos",
    }

    if reference_markers.issubset(cols):
        return "reference"

    observed_markers = {"Lof_varN", "syn_varN", "missense_positions"}
    if observed_markers.issubset(cols) and ("ENS_ID" in cols or "#ENS_ID" in cols):
        return "observed"

    raise ValueError(
        f"Could not detect whether {sumstats_file} is observed or reference-annotated. "
        "Use --sumstats_stage observed or --sumstats_stage reference."
    )


def prepare_reference_annotated_inputs(
    input_files,
    out_prefix: str,
    annotation_file: str | None = None,
    by_gene_MR_file: str | None = None,
    ENS2GeneID_file: str | None = None,
    MR_dist_file: str | None = None,
    sumstats_stage: str = "auto",
    write_variant_tables: bool = True,
    verbose: bool = True,
) -> list[str]:
    """
    Convert all inputs to reference-annotated summary-stat files.
    """
    file_process_dict = split_input_files(input_files)
    reference_annotated_files = []

    if file_process_dict["vcfs"]:
        require_existing_file(annotation_file, "variant annotation/mutation-target file")
        require_reference_args(by_gene_MR_file, ENS2GeneID_file, MR_dist_file)

    for vcf_file in file_process_dict["vcfs"]:
        if verbose:
            print(f"Preprocessing VCF: {vcf_file}")

        observed_sumstats_file = preprocess_one_vcf(
            vcf_file=vcf_file,
            annotation_file=annotation_file,
            out_prefix=out_prefix,
            write_variant_table=write_variant_tables,
        )

        if verbose:
            print(f"Adding reference annotations: {observed_sumstats_file}")

        reference_file = add_reference_to_observed_sumstats(
            observed_sumstats_file=observed_sumstats_file,
            by_gene_MR_file=by_gene_MR_file,
            ENS2GeneID_file=ENS2GeneID_file,
            MR_dist_file=MR_dist_file,
            out_prefix=out_prefix,
        )
        reference_annotated_files.append(reference_file)

    for sumstats_file in file_process_dict["sumstats"]:
        stage = detect_sumstats_stage(sumstats_file) if sumstats_stage == "auto" else sumstats_stage

        if stage == "reference":
            if verbose:
                print(f"Using reference-annotated summary stats: {sumstats_file}")
            reference_annotated_files.append(sumstats_file)

        elif stage == "observed":
            require_reference_args(by_gene_MR_file, ENS2GeneID_file, MR_dist_file)

            if verbose:
                print(f"Adding reference annotations to observed summary stats: {sumstats_file}")

            reference_file = add_reference_to_observed_sumstats(
                observed_sumstats_file=sumstats_file,
                by_gene_MR_file=by_gene_MR_file,
                ENS2GeneID_file=ENS2GeneID_file,
                MR_dist_file=MR_dist_file,
                out_prefix=out_prefix,
            )
            reference_annotated_files.append(reference_file)

        else:
            raise ValueError("--sumstats_stage must be one of: auto, observed, reference.")

    if not reference_annotated_files:
        raise ValueError("No usable input files were supplied.")

    if verbose:
        print("Working with the following reference-annotated summary-stat files:")
        for f in reference_annotated_files:
            print(f"  {f}")

    return reference_annotated_files


def combine_two_values(val1, val2):
    """
    Sum two scalar observed-stat values while preserving original NA behavior.
    """
    if pd.isna(val1) and pd.isna(val2):
        return np.nan
    if pd.isna(val1):
        return val2
    if pd.isna(val2):
        return val1
    return val1 + val2


def concat_two_cells(val1, val2):
    """
    Concatenate two comma-separated cells while preserving NA behavior.
    """
    if pd.isna(val1) and pd.isna(val2):
        return np.nan
    if pd.isna(val1):
        return val1 if False else val2
    if pd.isna(val2):
        return val1

    return ",".join(str(val1).split(",") + str(val2).split(","))


def join2_sumstat_dfs(df1, df2):
    """
    Join two reference-annotated summary-stat DataFrames.
    """
    validate_columns(df1.columns, {"ENS_ID"}, "first summary-stat DataFrame")
    validate_columns(df2.columns, {"ENS_ID"}, "second summary-stat DataFrame")

    df1_sorted = df1.sort_values(by="ENS_ID").reset_index(drop=True)
    df2_sorted = df2.sort_values(by="ENS_ID").reset_index(drop=True)

    if len(df1_sorted) != len(df2_sorted):
        raise ValueError(
            "Cannot join summary-stat files with different numbers of genes: "
            f"{len(df1_sorted)} vs {len(df2_sorted)}"
        )

    if not df1_sorted["ENS_ID"].equals(df2_sorted["ENS_ID"]):
        raise ValueError("Cannot join summary-stat files: ENS_ID columns do not match after sorting.")

    constant_columns = [
        "ENS_ID",
        "Gene_Name",
        "AM_lambda",
        "REVEL_lambda",
        "PAI_lambda",
        "Lof_lambda",
        "syn_lambda",
        "missense_lambda",
    ]
    sum_columns = ["AM_y", "REVEL_y", "PAI_y", "Lof_varN", "syn_varN"]
    concat_columns = ["missense_positions", "missense_MTs", "missense_MTs_by_pos"]

    validate_columns(
        df1_sorted.columns,
        constant_columns + sum_columns + concat_columns,
        "first summary-stat DataFrame",
    )
    validate_columns(
        df2_sorted.columns,
        constant_columns + sum_columns + concat_columns,
        "second summary-stat DataFrame",
    )

    sum_df = pd.DataFrame({})

    for col in constant_columns:
        sum_df[col] = df1_sorted[col]

    for col in sum_columns:
        sum_df[col] = [
            combine_two_values(val1, val2)
            for val1, val2 in zip(df1_sorted[col], df2_sorted[col])
        ]

    for col in concat_columns:
        sum_df[col] = [
            concat_two_cells(val1, val2)
            for val1, val2 in zip(df1_sorted[col], df2_sorted[col])
        ]

    return sum_df


def join_sumstat_files(sumstat_filelist):
    """
    Join multiple reference-annotated summary-stat files into one DataFrame.
    """
    if not sumstat_filelist:
        raise ValueError("sumstat_filelist is empty.")

    current_df = read_tsv(sumstat_filelist[0])

    for sumstat_file in sumstat_filelist[1:]:
        next_df = read_tsv(sumstat_file)
        current_df = join2_sumstat_dfs(current_df, next_df)

    return current_df


def import_optional_combine_module(module_name: str = "combine"):
    """
    Import combine.py if available, or return None.
    """
    try:
        if __package__:
            try:
                return importlib.import_module(f".{module_name}", package=__package__)
            except ImportError:
                pass
        return importlib.import_module(module_name)
    except ImportError:
        return None


def add_combined_pvalues(
    df: pd.DataFrame,
    combine_module_name: str = "combine",
    allow_legacy_pval_comb_lib: bool = True,
) -> pd.DataFrame:
    """
    Add ASCEND combined P-values.

    Preferred behavior:
        use combine.add_combined_pvalues(df), once combine.py exists.

    Temporary compatibility behavior:
        if combine.py is not available but pval_comb_lib.py is available, use
        its CCT/FCT functions to reproduce the original columns.
    """
    combine_module = import_optional_combine_module(combine_module_name)

    if combine_module is not None and hasattr(combine_module, "add_combined_pvalues"):
        return combine_module.add_combined_pvalues(df)

    if combine_module is None or not (hasattr(combine_module, "CCT") and hasattr(combine_module, "FCT")):
        if allow_legacy_pval_comb_lib:
            combine_module = import_optional_combine_module("pval_comb_lib")

    if combine_module is None or not (hasattr(combine_module, "CCT") and hasattr(combine_module, "FCT")):
        raise ImportError(
            "Could not combine P-values because neither combine.py with "
            "add_combined_pvalues(...) nor a module with CCT(...) and FCT(...) was found. "
            "Run with --skip_combine or add combine.py/pval_comb_lib.py."
        )

    out_df = df.copy()

    mis_enrich_Pvals = []
    for score_ID in ["AM", "REVEL", "PAI"]:
        col = f"{score_ID}_MisEnrich_Pval"
        validate_columns(out_df.columns, {col}, "component P-value DataFrame")
        mis_enrich_Pvals.append(out_df[col].values)

    mis_enrich_Pvals_T = np.array(mis_enrich_Pvals).T
    out_df["CCT_mis_enrich_Pval"] = combine_module.CCT(mis_enrich_Pvals_T)

    comb_Pvals = []
    for col in ["LoF_enrich_Pval", "CCT_mis_enrich_Pval", "Missense_Clust_Pval"]:
        validate_columns(out_df.columns, {col}, "component P-value DataFrame")
        comb_Pvals.append(out_df[col].values)

    out_df["FCT_all_Pval"] = combine_module.FCT(np.array(comb_Pvals).T)
    return out_df


def add_fdr_if_requested(
    df: pd.DataFrame,
    gene_constraint_file: str | None = None,
    OMIM_genes_file: str | None = None,
    censor_genes_file: str | None = None,
    skip_fdr: bool = False,
    normalize_weights: bool = False,
    add_reject_flags: bool = False,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Add wFDR columns unless skipped.
    """
    if skip_fdr:
        return df

    if fdr is None:
        raise ImportError("fdr.py could not be imported. Run with --skip_fdr or add fdr.py.")

    out_df = df.copy()

    if gene_constraint_file is not None:
        require_existing_file(gene_constraint_file, "gene wFDR weight file")
        out_df = fdr.add_gene_constraint_to_df(
            out_df,
            gene_constraint_file=gene_constraint_file,
            OMIM_genes_file=OMIM_genes_file,
        )
    else:
        validate_columns(out_df.columns, {"wFDR_weight"}, "results DataFrame")

    censor_col = None
    if censor_genes_file is not None:
        require_existing_file(censor_genes_file, "censor genes file")
        censor_col = "censor_flag"
        out_df = fdr.add_censor_flag_to_df(
            out_df,
            censor_genes_file=censor_genes_file,
            output_censor_col=censor_col,
        )

    out_df = fdr.add_fdr_columns(
        out_df,
        p_col="FCT_all_Pval",
        weight_col="wFDR_weight",
        normalize_weights=normalize_weights,
        censor_col=censor_col,
        alpha=alpha,
    )

    return out_df


def run_pipeline(
    input_files,
    out_prefix: str,
    annotation_file: str | None = None,
    by_gene_MR_file: str | None = None,
    ENS2GeneID_file: str | None = None,
    MR_dist_file: str | None = None,
    sumstats_stage: str = "auto",
    write_variant_tables: bool = True,
    syn_N: float | None = None,
    na_for0s: bool = False,
    clust_method: str = "gamma",
    clip_n_by: int = 2,
    clip_clust_pvals_by: float = 1e-10,
    suppress_clust: bool = False,
    Uprod_file: str | None = None,
    skip_combine: bool = False,
    combine_module: str = "combine",
    skip_fdr: bool = False,
    gene_constraint_file: str | None = None,
    OMIM_genes_file: str | None = None,
    censor_genes_file: str | None = None,
    normalize_weights: bool = False,
    add_reject_flags: bool = False,
    alpha: float = 0.05,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Run the end-to-end ASCEND pipeline.
    """
    reference_files = prepare_reference_annotated_inputs(
        input_files=input_files,
        out_prefix=out_prefix,
        annotation_file=annotation_file,
        by_gene_MR_file=by_gene_MR_file,
        ENS2GeneID_file=ENS2GeneID_file,
        MR_dist_file=MR_dist_file,
        sumstats_stage=sumstats_stage,
        write_variant_tables=write_variant_tables,
        verbose=verbose,
    )

    if verbose:
        print("Joining summary-stat files...")

    sumstat_df = join_sumstat_files(reference_files)

    joined_sumstats_file = f"{out_prefix}_sumstats.txt"
    write_tsv(sumstat_df, joined_sumstats_file)

    if verbose:
        print(f"Wrote joined summary stats: {joined_sumstats_file}")
        print(f"Total synonymous variant N: {sumstat_df['syn_varN'].sum()}")

    if clust_method == "permut":
        require_existing_file(Uprod_file, "Uprod clustering distribution file")

    if verbose:
        print("Calculating component association statistics...")

    stats_df = statistics.add_association_statistics(
        sumstat_df,
        syn_N=syn_N,
        na_for0s=na_for0s,
        clust_method=clust_method,
        clip_n_by=clip_n_by,
        clip_clust_pvals_by=clip_clust_pvals_by,
        suppress_clust=suppress_clust,
        Uprod_file=Uprod_file,
    )

    component_pvals_file = f"{out_prefix}_component_pvals.tsv"
    write_tsv(stats_df, component_pvals_file)

    if verbose:
        print(f"Wrote component P-values: {component_pvals_file}")

    results_df = stats_df

    if not skip_combine:
        if verbose:
            print("Combining component P-values...")
        results_df = add_combined_pvalues(results_df, combine_module_name=combine_module)

    if not skip_fdr:
        if verbose:
            print("Adding weighted/censored FDR columns...")
        results_df = add_fdr_if_requested(
            results_df,
            gene_constraint_file=gene_constraint_file,
            OMIM_genes_file=OMIM_genes_file,
            censor_genes_file=censor_genes_file,
            skip_fdr=skip_fdr,
            normalize_weights=normalize_weights,
            add_reject_flags=add_reject_flags,
            alpha=alpha,
        )

    if "FCT_all_Qval" in results_df.columns:
        results_df = results_df.sort_values(by="FCT_all_Qval", na_position="last").reset_index(drop=True)
    elif "FCT_all_Pval" in results_df.columns:
        results_df = results_df.sort_values(by="FCT_all_Pval", na_position="last").reset_index(drop=True)

    final_results_file = f"{out_prefix}_sumstats_results.txt"
    write_tsv(results_df, final_results_file)

    if verbose:
        print(f"Wrote final results: {final_results_file}")

    return results_df


if __name__ == "__main__":

    EXAMPLES = """
Basic use
---------
  python pipeline.py --i cohort.vcf.gz --o ASCEND

Multiple cohorts
----------------
  python pipeline.py --i cohort1.vcf.gz,cohort2.vcf.gz --o ASCEND

Use a different packaged data directory
---------------------------------------
  python pipeline.py --i cohort.vcf.gz --o ASCEND --data_dir /path/to/data

Show all advanced options
-------------------------
  python pipeline.py --advanced-help
"""

    def add_argument(arg_container, *args, advanced: bool = False, **kwargs):
        """
        Add an argparse option.

        Advanced options are still accepted by the CLI, but are hidden from the
        default help menu. They appear only with --advanced-help.
        """
        show_advanced = getattr(arg_container, "_ascend_show_advanced_options", False)
        if advanced and not show_advanced:
            kwargs["help"] = argparse.SUPPRESS
        arg_container.add_argument(*args, **kwargs)

    class ASCENDArgumentParser(argparse.ArgumentParser):
        def __init__(self, *args, show_advanced_options: bool = False, **kwargs):
            self.show_advanced_options = show_advanced_options
            super().__init__(*args, **kwargs)

    def build_parser(show_advanced_options: bool = False):
        parser = ASCENDArgumentParser(
            description=(
                "Run ASCEND on one or more VCF files. In the standard setup, "
                "all reference files are read from ../data relative to pipeline.py, "
                "so only --i and --o are required."
            ),
            epilog=EXAMPLES,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            show_advanced_options=show_advanced_options,
        )

        required = parser.add_argument_group("required arguments")
        required.add_argument(
            "--i",
            required=True,
            help="Input VCF, or comma-separated list of VCFs. Summary-stat files are also accepted.",
        )
        required.add_argument(
            "--o",
            required=True,
            help="Output prefix.",
        )

        common = parser.add_argument_group("common optional arguments")
        common.add_argument(
            "--data_dir",
            default=None,
            help="Directory containing ASCEND reference files. Default: ../data relative to pipeline.py.",
        )
        common.add_argument(
            "--sumstats_stage",
            default="auto",
            choices=["auto", "observed", "reference"],
            help="How to treat non-VCF inputs. Default: auto.",
        )
        common.add_argument(
            "--quiet",
            action="store_true",
            help="Reduce progress printing.",
        )
        common.add_argument(
            "--advanced-help",
            action="store_true",
            help="Show all configuration and method-development options.",
        )

        references = parser.add_argument_group("reference file overrides")
        references._ascend_show_advanced_options = show_advanced_options
        add_argument(
            references,
            "--annotation",
            default=None,
            advanced=True,
            help="Override variant annotation/mutation-target table.",
        )
        add_argument(
            references,
            "--by_gene_MR",
            default=None,
            advanced=True,
            help="Override by-gene mutation-rate reference file.",
        )
        add_argument(
            references,
            "--ENS2GeneID",
            default=None,
            advanced=True,
            help="Override ENS_ID-to-gene-name mapping file.",
        )
        add_argument(
            references,
            "--MR_dist",
            default=None,
            advanced=True,
            help="Override per-position MR-distance reference file.",
        )
        add_argument(
            references,
            "--Uprod_file",
            default=None,
            advanced=True,
            help="Override permutation clustering distribution file.",
        )

        output = parser.add_argument_group("output controls")
        output._ascend_show_advanced_options = show_advanced_options
        add_argument(
            output,
            "--no_variant_tables",
            action="store_true",
            advanced=True,
            help="Do not write observed variant tables for VCF inputs.",
        )

        stats = parser.add_argument_group("statistical options")
        stats._ascend_show_advanced_options = show_advanced_options
        add_argument(
            stats,
            "--syn_N",
            type=float,
            default=None,
            advanced=True,
            help="Total synonymous count for mutation-rate scaling. Default: sum(syn_varN) after joining.",
        )
        add_argument(
            stats,
            "--r",
            action="store_true",
            advanced=True,
            help="Convert LoF P-values equal to 1 to NA before P-value combination.",
        )
        add_argument(
            stats,
            "--clust_method",
            default="gamma",
            choices=["gamma", "chi2", "permut"],
            advanced=True,
            help="Method for clustering P-values. Default: gamma.",
        )
        add_argument(
            stats,
            "--m",
            type=int,
            default=2,
            advanced=True,
            help="Minimum number of missense variants for clustering. Default: 2.",
        )
        add_argument(
            stats,
            "--c",
            type=float,
            default=1e-10,
            advanced=True,
            help="Lower clipping threshold for clustering P-values. Default: 1e-10.",
        )
        add_argument(
            stats,
            "--suppress_clust",
            action="store_true",
            advanced=True,
            help="Do not calculate the clustering statistic.",
        )

        combine_group = parser.add_argument_group("P-value combination options")
        combine_group._ascend_show_advanced_options = show_advanced_options
        add_argument(
            combine_group,
            "--skip_combine",
            action="store_true",
            advanced=True,
            help="Stop after component P-values and do not add CCT/FCT combined P-values.",
        )
        add_argument(
            combine_group,
            "--combine_module",
            default="combine",
            advanced=True,
            help="Module used for P-value combination. Default: combine.",
        )

        fdr_group = parser.add_argument_group("FDR options")
        fdr_group._ascend_show_advanced_options = show_advanced_options
        add_argument(
            fdr_group,
            "--skip_fdr",
            action="store_true",
            advanced=True,
            help="Do not add weighted/censored FDR columns.",
        )
        add_argument(
            fdr_group,
            "--gene_constraint_file",
            default=None,
            advanced=True,
            help="Override gene wFDR weight file.",
        )
        add_argument(
            fdr_group,
            "--OMIM_genes_file",
            default=None,
            advanced=True,
            help="Override one-column ENS_ID file used to add OMIM_flag.",
        )
        add_argument(
            fdr_group,
            "--censor_genes_file",
            default=None,
            advanced=True,
            help="One-column ENS_ID file for censored BH.",
        )
        add_argument(
            fdr_group,
            "--normalize_weights",
            action="store_true",
            advanced=True,
            help="Normalize weights to mean 1 over tested genes.",
        )
        add_argument(
            fdr_group,
            "--add_reject_flags",
            action="store_true",
            advanced=True,
            help="Also add Boolean rejection columns at --alpha.",
        )
        add_argument(
            fdr_group,
            "--alpha",
            type=float,
            default=0.05,
            advanced=True,
            help="FDR threshold used only with --add_reject_flags. Default: 0.05.",
        )

        return parser

    if "--advanced-help" in sys.argv:
        advanced_parser = build_parser(show_advanced_options=True)
        advanced_parser.print_help()
        sys.exit(0)

    parser = build_parser(show_advanced_options=False)
    args = parser.parse_args()

    default_paths = resolve_data_paths(args.data_dir)

    annotation_file = use_default_if_none(args.annotation, default_paths["annotation_file"])
    by_gene_MR_file = use_default_if_none(args.by_gene_MR, default_paths["by_gene_MR_file"])
    ENS2GeneID_file = use_default_if_none(args.ENS2GeneID, default_paths["ENS2GeneID_file"])
    MR_dist_file = use_default_if_none(args.MR_dist, default_paths["MR_dist_file"])
    Uprod_file = use_default_if_none(args.Uprod_file, default_paths["Uprod_file"])

    gene_constraint_file = use_default_if_none(
        args.gene_constraint_file,
        default_paths["gene_constraint_file"],
    )
    OMIM_genes_file = use_default_if_none(
        args.OMIM_genes_file,
        default_paths["OMIM_genes_file"],
    )

    input_files = [x.strip() for x in args.i.split(",") if x.strip()]

    if not args.quiet:
        print(f"ASCEND data directory: {Path(default_paths['annotation_file']).parent}")

    run_pipeline(
        input_files=input_files,
        out_prefix=args.o,
        annotation_file=annotation_file,
        by_gene_MR_file=by_gene_MR_file,
        ENS2GeneID_file=ENS2GeneID_file,
        MR_dist_file=MR_dist_file,
        sumstats_stage=args.sumstats_stage,
        write_variant_tables=not args.no_variant_tables,
        syn_N=args.syn_N,
        na_for0s=args.r,
        clust_method=args.clust_method,
        clip_n_by=args.m,
        clip_clust_pvals_by=args.c,
        suppress_clust=args.suppress_clust,
        Uprod_file=Uprod_file,
        skip_combine=args.skip_combine,
        combine_module=args.combine_module,
        skip_fdr=args.skip_fdr,
        gene_constraint_file=gene_constraint_file,
        OMIM_genes_file=OMIM_genes_file,
        censor_genes_file=args.censor_genes_file,
        normalize_weights=args.normalize_weights,
        add_reject_flags=args.add_reject_flags,
        alpha=args.alpha,
        verbose=not args.quiet,
    )

