"""
combine.py

P-value combination utilities for ASCEND.

This module implements the two P-value combination steps used by the current
ASCEND:

1. CCT / ACAT-style Cauchy combination for correlated missense-score P-values.
2. Fisher's combination test for combining broader evidence classes.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import cauchy, chi2


def CCT_ind(p_values, weights_vec = None):
    # Check for NaN values
    p_values = np.array([p for p in p_values if not np.isnan(p)])
    if 0 in p_values:
        return 0
    
    # Compute/normalize the weights
    if weights_vec == None:
        l = len(p_values)
        weights_vec = np.array([1.0/l for i in range(l)])
    
    weights_vec = np.array(weights_vec)/np.sum(weights_vec)
    
    # Apply tangent transformation to P-values
    transformed_p_values = np.tan(np.pi * (0.5 - p_values))

    # Calculate the test statistic
    test_statistic = np.sum(weights_vec*transformed_p_values)

    # Calculate the combined P-value
    combined_p_value = cauchy.sf(test_statistic)

    return combined_p_value


def CCT(pvals, weights_vec = None):
    # pvals : list of arrays of Ps to be combined
    pvals = np.array(pvals)
    if pvals.ndim == 1:
        return CCT_ind(pvals, weights_vec)
    
    return np.array([CCT_ind(pvals[i], weights_vec) for i in range(len(pvals))])


def FCT_ind(p_values):
    # Check for NaN values
    p_values = np.array([p for p in p_values if not np.isnan(p)])
    if 0 in p_values:
        return 0
    
    L = len(p_values)
    return chi2.sf(np.sum(-2*np.log(p_values)), 2*L)


def FCT(pvals):
    pvals = np.array(pvals)
    if pvals.ndim == 1:
        return FCT_ind(pvals)

    return np.array([FCT_ind(pvals[i]) for i in range(len(pvals))])


def combine_columns(
    df: pd.DataFrame,
    columns: list[str],
    method: str,
    output_column: str,
    weights_vec=None,
) -> pd.DataFrame:
    """
    Add a combined-P-value column to a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input gene-level table.
    columns : list[str]
        P-value columns to combine row-wise.
    method : {"CCT", "FCT"}
        Combination method.
    output_column : str
        Name of the new combined-P-value column.
    weights_vec : array-like or None
        Optional weights for CCT. Ignored for FCT.

    Returns
    -------
    pd.DataFrame
        Copy of df with output_column added.
    """
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"Input DataFrame is missing columns: {missing}")

    method = method.upper()
    pval_matrix = df[columns].to_numpy(dtype=float)

    out_df = df.copy()

    if method == "CCT":
        out_df[output_column] = CCT(pval_matrix, weights_vec=weights_vec)
    elif method == "FCT":
        out_df[output_column] = FCT(pval_matrix)
    else:
        raise ValueError("method must be one of: 'CCT', 'FCT'")

    return out_df


def add_missense_cct(
    df: pd.DataFrame,
    score_ids: list[str] | None = None,
    output_column: str = "CCT_mis_enrich_Pval",
    weights_vec=None,
) -> pd.DataFrame:
    """
    Combine missense enrichment P-values across score models using CCT.

    By default, combines:
        AM_MisEnrich_Pval
        REVEL_MisEnrich_Pval
        PAI_MisEnrich_Pval
    """
    if score_ids is None:
        score_ids = ["AM", "REVEL", "PAI"]

    columns = [f"{score_id}_MisEnrich_Pval" for score_id in score_ids]
    return combine_columns(
        df=df,
        columns=columns,
        method="CCT",
        output_column=output_column,
        weights_vec=weights_vec,
    )


def add_fisher_all(
    df: pd.DataFrame,
    columns: list[str] | None = None,
    output_column: str = "FCT_all_Pval",
) -> pd.DataFrame:
    """
    Combine the main ASCEND evidence classes using Fisher's combination test.

    By default, combines:
        LoF_enrich_Pval
        CCT_mis_enrich_Pval
        Missense_Clust_Pval

    NaNs are omitted row by row, matching the current ASCEND missing-data
    behavior.
    """
    if columns is None:
        columns = [
            "LoF_enrich_Pval",
            "CCT_mis_enrich_Pval",
            "Missense_Clust_Pval",
        ]

    return combine_columns(
        df=df,
        columns=columns,
        method="FCT",
        output_column=output_column,
    )
