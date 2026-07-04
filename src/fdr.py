"""
fdr.py

Weighted and censored Benjamini-Hochberg utilities for ASCEND.


ASCEND convention
-----------------
For a gene g with P-value p_g and external nonnegative weight w_g, the weighted
P-value used for weighted BH is

    Q_g = p_g / w_g.

The weights are assumed to be independent of the de novo P-values under the
null. In the manuscript, weights are normalized to have mean 1 across tested
genes. This module can normalize weights if requested, but by default it uses
the weights as supplied, matching the original analysis script.

Censoring convention
--------------------
Censored genes are removed from the BH procedure entirely. They keep their
weighted P-values for inspection, but their BH q-values are set to NA and their
rejection flags are False.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd


def read_table(path_str: str, **kwargs) -> pd.DataFrame:
    """Read a TSV/TSV.GZ file with compression inferred from the filename."""
    return pd.read_table(path_str, sep="\t", compression="infer", **kwargs)


def validate_columns(columns, required, source_name: str = "") -> None:
    """Raise a readable error if required columns are missing."""
    missing = set(required) - set(columns)
    if missing:
        label = f" in {source_name}" if source_name else ""
        raise ValueError(f"Missing required columns{label}: {sorted(missing)}")


def validate_pvalues(p_values, allow_nan: bool = True) -> np.ndarray:
    """
    Convert input to a float array and check that finite P-values are in [0, 1].
    """
    p = np.asarray(p_values, dtype=float)
    finite = np.isfinite(p)

    if not allow_nan and np.any(~finite):
        raise ValueError("P-values contain NaN or infinite values.")

    bad = finite & ((p < 0.0) | (p > 1.0))
    if np.any(bad):
        example = p[np.where(bad)[0][0]]
        raise ValueError(f"P-values must be in [0, 1]. Found value: {example}")

    return p


def validate_weights(weights, allow_nan: bool = False) -> np.ndarray:
    """
    Convert input to a float array and check that finite weights are positive.
    """
    w = np.asarray(weights, dtype=float)
    finite = np.isfinite(w)

    if not allow_nan and np.any(~finite):
        raise ValueError("Weights contain NaN or infinite values.")

    bad = finite & (w <= 0.0)
    if np.any(bad):
        example = w[np.where(bad)[0][0]]
        raise ValueError(f"Weights must be positive. Found value: {example}")

    return w


def normalize_weights_to_mean_one(weights, mask=None) -> np.ndarray:
    """
    Normalize weights so that their mean is 1 over the selected genes.

    Parameters
    ----------
    weights : array-like
        Positive weights.
    mask : array-like of bool, optional
        If supplied, normalization is performed over mask == True only. Values
        outside the mask are still scaled by the same normalization constant.

    Returns
    -------
    np.ndarray
        Scaled weights.
    """
    w = validate_weights(weights, allow_nan=False)

    if mask is None:
        norm_mask = np.ones(len(w), dtype=bool)
    else:
        norm_mask = np.asarray(mask, dtype=bool)
        if len(norm_mask) != len(w):
            raise ValueError("mask and weights have different lengths.")

    if not np.any(norm_mask):
        raise ValueError("Cannot normalize weights: mask selects zero genes.")

    mean_w = np.mean(w[norm_mask])
    if mean_w <= 0.0 or not np.isfinite(mean_w):
        raise ValueError(f"Cannot normalize weights: invalid mean weight {mean_w}")

    return w / mean_w


def weighted_pvalues(p_values, weights, normalize_weights: bool = False, normalize_mask=None):
    """
    Calculate weighted P-values p_g / w_g.

    Missing P-values remain NaN. This function does not impute missing values
    to 1, because downstream ASCEND combination code treats missing components
    as omitted rather than neutralized.
    """
    p = validate_pvalues(p_values, allow_nan=True)
    w = validate_weights(weights, allow_nan=False)

    if len(p) != len(w):
        raise ValueError("p_values and weights have different lengths.")

    if normalize_weights:
        w = normalize_weights_to_mean_one(w, mask=normalize_mask)

    out = np.full(len(p), np.nan, dtype=float)
    valid = np.isfinite(p)
    out[valid] = p[valid] / w[valid]
    return out


def bh_adjust(p_values, mask=None) -> np.ndarray:
    """
    Benjamini-Hochberg adjusted q-values.

    Parameters
    ----------
    p_values : array-like
        P-values to adjust. NaN values are ignored and remain NaN.
    mask : array-like of bool, optional
        If supplied, only entries with mask == True are included in the BH
        procedure. Entries outside the mask remain NaN.

    Returns
    -------
    np.ndarray
        BH q-values, aligned to the original input order.
    """
    p = validate_pvalues(p_values, allow_nan=True)

    if mask is None:
        test_mask = np.ones(len(p), dtype=bool)
    else:
        test_mask = np.asarray(mask, dtype=bool)
        if len(test_mask) != len(p):
            raise ValueError("mask and p_values have different lengths.")

    valid = test_mask & np.isfinite(p)
    out = np.full(len(p), np.nan, dtype=float)

    if not np.any(valid):
        return out

    valid_idx = np.where(valid)[0]
    p_valid = p[valid]

    order = np.argsort(p_valid)
    p_sorted = p_valid[order]
    m = len(p_sorted)
    ranks = np.arange(1, m + 1, dtype=float)

    q_sorted = p_sorted * m / ranks
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.minimum(q_sorted, 1.0)

    out[valid_idx[order]] = q_sorted
    return out


def bh_reject(p_values, alpha: float = 0.05, mask=None) -> np.ndarray:
    """
    Boolean rejection flags from BH-adjusted q-values.
    """
    if alpha < 0.0 or alpha > 1.0:
        raise ValueError("alpha must be in [0, 1].")

    q = bh_adjust(p_values, mask=mask)
    return np.isfinite(q) & (q <= alpha)


def weighted_bh_adjust(
    p_values,
    weights,
    normalize_weights: bool = False,
    mask=None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Weighted BH q-values.

    Returns
    -------
    weighted_p : np.ndarray
        p / weight, aligned to input order. May exceed 1.
    q_values : np.ndarray
        BH q-values obtained by applying BH to weighted_p over mask == True.
    """
    if mask is None:
        norm_mask = None
    else:
        norm_mask = np.asarray(mask, dtype=bool)

    weighted_p = weighted_pvalues(
        p_values,
        weights,
        normalize_weights=normalize_weights,
        normalize_mask=norm_mask,
    )

    # BH expects values in [0, 1], whereas weighted P-values can exceed 1.
    # Clipping above 1 does not affect rejections and gives standard q-values.
    weighted_p_for_bh = np.minimum(weighted_p, 1.0)
    q_values = bh_adjust(weighted_p_for_bh, mask=mask)

    return weighted_p, q_values


def censored_weighted_bh_adjust(
    p_values,
    weights,
    censored,
    normalize_weights: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Weighted BH after removing censored genes from the tested set.

    Parameters
    ----------
    censored : array-like of bool
        True means the gene is censored and is not included in the BH procedure.
    """
    censored = np.asarray(censored, dtype=bool)

    if len(censored) != len(np.asarray(p_values)):
        raise ValueError("censored and p_values have different lengths.")

    uncensored_mask = ~censored

    weighted_p, q_values = weighted_bh_adjust(
        p_values,
        weights,
        normalize_weights=normalize_weights,
        mask=uncensored_mask,
    )
    return weighted_p, q_values


def add_gene_constraint_to_df(
    merged_df: pd.DataFrame,
    gene_constraint_file: str,
    OMIM_genes_file: str | None = None,
    weight_file_column: str = "50_quant_weight",
    ens_col: str = "ENS_ID",
    output_weight_col: str = "wFDR_weight",
    output_omim_col: str = "OMIM_flag",
    default_weight: float = 1.0,
) -> pd.DataFrame:
    """
    Add ASCEND gene-level wFDR weights to a result DataFrame.

    This mirrors the original analysis-script behavior:
    - read weights from gene_constraint_file;
    - map by ENS_ID;
    - assign default_weight to genes missing from the weight file;
    - optionally add OMIM_flag from a one-column ENS_ID file.
    """
    validate_columns(merged_df.columns, {ens_col}, "analysis result DataFrame")

    gene_constraint_df = pd.read_csv(gene_constraint_file, sep="\t")
    validate_columns(
        gene_constraint_df.columns,
        {ens_col, weight_file_column},
        "gene constraint weight file",
    )

    weight_dict = {
        gene_constraint_df[ens_col].iloc[i]: gene_constraint_df[weight_file_column].iloc[i]
        for i in range(len(gene_constraint_df))
    }

    out_df = merged_df.copy()
    out_df[output_weight_col] = [
        weight_dict.get(out_df[ens_col].iloc[i], default_weight)
        for i in range(len(out_df))
    ]

    # Validate after assigning default weights.
    validate_weights(out_df[output_weight_col].values, allow_nan=False)

    if OMIM_genes_file is not None:
        omim_ens_ids = pd.read_table(OMIM_genes_file, header=None)[0].tolist()
        out_df[output_omim_col] = out_df[ens_col].isin(omim_ens_ids).astype(int)

    return out_df


def add_censor_flag_to_df(
    merged_df: pd.DataFrame,
    censor_genes_file: str,
    ens_col: str = "ENS_ID",
    output_censor_col: str = "censor_flag",
) -> pd.DataFrame:
    """
    Add a censoring flag from a one-column ENS_ID file.

    The censoring file should contain genes to remove from the censored BH
    procedure, one ENS_ID per line. A header is not required.
    """
    validate_columns(merged_df.columns, {ens_col}, "analysis result DataFrame")

    censor_ens_ids = pd.read_table(censor_genes_file, header=None)[0].tolist()

    out_df = merged_df.copy()
    out_df[output_censor_col] = out_df[ens_col].isin(censor_ens_ids).astype(int)
    return out_df


def add_fdr_columns(
    merged_df: pd.DataFrame,
    p_col: str = "FCT_all_Pval",
    weight_col: str = "wFDR_weight",
    alpha: float = 0.05,
    normalize_weights: bool = False,
    censor_col: str | None = None,
    legacy_weighted_p_col: str = "FCT_all_Qval",
    weighted_bh_q_col: str = "FCT_all_wBH_Qval",
    weighted_bh_reject_col: str = "FCT_all_wBH_reject",
    censored_weighted_bh_q_col: str = "FCT_all_cens_wBH_Qval",
    censored_weighted_bh_reject_col: str = "FCT_all_cens_wBH_reject",
) -> pd.DataFrame:
    """
    Add weighted and optionally censored weighted BH columns.

    The default `legacy_weighted_p_col` is intentionally named FCT_all_Qval to
    preserve the original analysis-script convention:

        FCT_all_Qval = FCT_all_Pval / wFDR_weight

    Actual BH-adjusted q-values are written to FCT_all_wBH_Qval and, if
    censor_col is supplied, FCT_all_cens_wBH_Qval.
    """
    validate_columns(merged_df.columns, {p_col, weight_col}, "analysis result DataFrame")

    out_df = merged_df.copy()

    p_values = out_df[p_col].values
    weights = out_df[weight_col].values

    weighted_p, weighted_q = weighted_bh_adjust(
        p_values,
        weights,
        normalize_weights=normalize_weights,
        mask=None,
    )

    out_df[legacy_weighted_p_col] = weighted_p
    out_df[weighted_bh_q_col] = weighted_q
    out_df[weighted_bh_reject_col] = np.isfinite(weighted_q) & (weighted_q <= alpha)

    if censor_col is not None:
        validate_columns(out_df.columns, {censor_col}, "analysis result DataFrame")
        censored = out_df[censor_col].astype(bool).values

        _, censored_q = censored_weighted_bh_adjust(
            p_values,
            weights,
            censored=censored,
            normalize_weights=normalize_weights,
        )

        out_df[censored_weighted_bh_q_col] = censored_q
        out_df[censored_weighted_bh_reject_col] = (
            np.isfinite(censored_q) & (censored_q <= alpha)
        )

    return out_df


def apply_fdr_from_files(
    input_file: str,
    output_file: str,
    gene_constraint_file: str | None = None,
    OMIM_genes_file: str | None = None,
    censor_genes_file: str | None = None,
    p_col: str = "FCT_all_Pval",
    weight_file_column: str = "50_quant_weight",
    weight_col: str = "wFDR_weight",
    alpha: float = 0.05,
    normalize_weights: bool = False,
) -> pd.DataFrame:
    """
    File-level wrapper for adding wFDR/censored-wFDR columns.
    """
    df = pd.read_csv(input_file, sep="\t")

    if gene_constraint_file is not None:
        df = add_gene_constraint_to_df(
            df,
            gene_constraint_file=gene_constraint_file,
            OMIM_genes_file=OMIM_genes_file,
            weight_file_column=weight_file_column,
            output_weight_col=weight_col,
        )
    else:
        validate_columns(df.columns, {weight_col}, "input file")

    censor_col = None
    if censor_genes_file is not None:
        censor_col = "censor_flag"
        df = add_censor_flag_to_df(
            df,
            censor_genes_file=censor_genes_file,
            output_censor_col=censor_col,
        )

    df = add_fdr_columns(
        df,
        p_col=p_col,
        weight_col=weight_col,
        alpha=alpha,
        normalize_weights=normalize_weights,
        censor_col=censor_col,
    )

    # Keep the original sorting convention: sort by weighted P-value.
    df = df.sort_values(by="FCT_all_Qval", na_position="last").reset_index(drop=True)
    df.to_csv(output_file, sep="\t", index=False, na_rep="NA")
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add ASCEND weighted BH and censored weighted BH columns."
    )
    parser.add_argument("--i", required=True, help="Input ASCEND result TSV")
    parser.add_argument("--o", required=True, help="Output TSV")
    parser.add_argument("--p_col", default="FCT_all_Pval", help="Column containing gene-level P-values")
    parser.add_argument("--weight_col", default="wFDR_weight", help="Column containing gene weights")
    parser.add_argument(
        "--gene_constraint_file",
        default=None,
        help="Optional gene weight file. If supplied, weights are added before FDR calculation.",
    )
    parser.add_argument(
        "--weight_file_column",
        default="50_quant_weight",
        help="Weight column in --gene_constraint_file.",
    )
    parser.add_argument(
        "--OMIM_genes_file",
        default=None,
        help="Optional one-column ENS_ID file used to add OMIM_flag.",
    )
    parser.add_argument(
        "--censor_genes_file",
        default=None,
        help="Optional one-column ENS_ID file for censored BH.",
    )
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR level for rejection flags")
    parser.add_argument(
        "--normalize_weights",
        action="store_true",
        help="Normalize weights to mean 1 over tested genes before applying weighted BH.",
    )
    args = parser.parse_args()

    apply_fdr_from_files(
        input_file=args.i,
        output_file=args.o,
        gene_constraint_file=args.gene_constraint_file,
        OMIM_genes_file=args.OMIM_genes_file,
        censor_genes_file=args.censor_genes_file,
        p_col=args.p_col,
        weight_file_column=args.weight_file_column,
        weight_col=args.weight_col,
        alpha=args.alpha,
        normalize_weights=args.normalize_weights,
    )
