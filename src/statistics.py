"""
statistics.py

Association-statistic utilities for ASCEND.

This module calculates the gene-level statistical components used by ASCEND:

1. LoF enrichment P-values using a one-sided Poisson tail;
2. missense enrichment P-values using the Poisson-Irwin-Hall statistic;
3. missense clustering P-values using the restricted spacing statistic.

The function names and output column names are kept close to the original
analysis script to avoid changing downstream behavior.
"""

from __future__ import annotations
import argparse
import gzip
import numpy as np
import pandas as pd
from scipy.stats import chi2, gamma, irwinhall, poisson


def read_table(path_str: str, **kwargs) -> pd.DataFrame:
    """Read a TSV/TSV.GZ file with compression inferred from the filename."""
    return pd.read_table(path_str, sep="\t", compression="infer", **kwargs)


def open_text(path_str: str):
    """Open plain-text or gzipped files for reading."""
    if path_str.endswith(".gz"):
        return gzip.open(path_str, "rt")
    return open(path_str, "r")


def validate_columns(columns, required, source_name: str = "") -> None:
    """Raise a readable error if required columns are missing."""
    missing = set(required) - set(columns)
    if missing:
        label = f" in {source_name}" if source_name else ""
        raise ValueError(f"Missing required columns{label}: {sorted(missing)}")


def calc_LoF_enrich_P(merged_df, syn_N, na_for0s=False):
    """
    Calculate LoF enrichment P-values using a one-sided Poisson survival function.

        lambda_g = Lof_lambda_g * total_synonymous_count
        P_g = P{Pois(lambda_g) >= observed_LoF_count_g}

    Parameters
    ----------
    merged_df : pd.DataFrame
        Reference-annotated summary-stat DataFrame.
    syn_N : float
        Total observed synonymous count used for cohort-level scaling.
    na_for0s : bool
        If True, P-values >= 1 are converted to NaN
    """
    validate_columns(
        merged_df.columns,
        {"Lof_lambda", "Lof_varN"},
        "reference-annotated summary-stat DataFrame",
    )

    merged_df = merged_df.copy()

    LoF_lambda = merged_df["Lof_lambda"] * syn_N
    merged_df["LoF_enrich_Pval"] = poisson.sf(merged_df["Lof_varN"] - 1, LoF_lambda)

    if na_for0s:
        merged_df.loc[merged_df["LoF_enrich_Pval"] >= 1.0, "LoF_enrich_Pval"] = np.nan

    return merged_df


def PoissonIH_SF(y, lam):
    """
    Calculate the Poisson-Irwin-Hall survival function.

    The statistic is a mixture over the number of observed missense variants:
    for i >= 1, the conditional score-sum distribution is Irwin-Hall(i), and
    the mixing weights are Poisson(i | lambda).

    This function preserves the truncation rule from the original analysis
    script. Missing y/lambda values return NaN.
    """
    if pd.isna(lam) or pd.isna(y):
        return np.nan
    if lam == 0:
        return np.nan

    P_val = 0.0
    summation_range = max(
        20,
        int(lam + np.sqrt(lam) * 20),
        int(y + np.sqrt(y) * 10),
    )

    for i in range(1, summation_range):
        pois_p = poisson.pmf(i, lam)
        ih_sf = irwinhall.sf(y, i)
        P_val += pois_p * ih_sf

    if P_val > 1.0:
        print(f"Warning: Poisson-IH SF P-value > 1.0 : {P_val}", y, lam)
    if P_val <= 0.0:
        print(f"Warning: Poisson-IH SF P-value <= 0.0 : {P_val}", y, lam)

    return P_val


def calc_MisEnrich_P(merged_df, score_ID, synN):
    """
    Calculate missense enrichment P-values for one functionality score.

    Parameters
    ----------
    merged_df : pd.DataFrame
        Reference-annotated summary-stat DataFrame.
    score_ID : str
        One of "AM", "REVEL", "PAI" by default.
    synN : float
        Total observed synonymous count used for cohort-level scaling.

    Adds
    ----
    {score_ID}_MisEnrich_Pval
    """
    lambda_column = f"{score_ID}_lambda"
    y_column = f"{score_ID}_y"
    output_column = f"{score_ID}_MisEnrich_Pval"

    validate_columns(
        merged_df.columns,
        {lambda_column, y_column},
        "reference-annotated summary-stat DataFrame",
    )

    merged_df = merged_df.copy()

    P_vals = []
    for _, row in merged_df.iterrows():
        lam = row[lambda_column] * synN
        y = row[y_column]

        if pd.isna(lam) or pd.isna(y):
            P_vals.append(np.nan)
        else:
            P_vals.append(PoissonIH_SF(y, lam))

    merged_df[output_column] = P_vals
    return merged_df


class SimplexVolumeDist:
    """
    Utilities for the ASCEND missense clustering statistic.

    The default clustering method is the shifted Gamma approximation for the
    restricted internal-spacing statistic S'. The optional permutation-table
    method is supported if a Uprod distribution file is supplied.
    """

    def __init__(self, Uprod_file: str | None = None):
        self.Uprod_file = Uprod_file
        self.D_x = {}
        self.D_CDF = {}
        self.L = 0

        if Uprod_file is not None:
            self.load_Uprod_file(Uprod_file)

    def load_Uprod_file(self, Uprod_file: str) -> None:
        """
        Load precomputed permutation/empirical CDF tables for the clustering statistic.

        Expected format:
            n_points    comma_separated_x_values    comma_separated_cdf_values
        """
        D_x = {}
        D_CDF = {}

        with open_text(Uprod_file) as infile:
            for dist_str in infile:
                dist_str = dist_str.strip().split("\t")
                n = int(dist_str[0])

                x_arr = np.array(list(map(float, dist_str[1].split(","))))
                cdf_arr = np.array(list(map(float, dist_str[2].split(","))))

                D_x[n] = x_arr
                D_CDF[n] = cdf_arr

        self.Uprod_file = Uprod_file
        self.D_x = D_x
        self.D_CDF = D_CDF
        self.L = len(cdf_arr) if D_CDF else 0

    def harmonic_numbers(self, n, power_fac=1):
        """
        Harmonic numbers of order n and power_fac:
            H_n^(power_fac) = sum_{k=1}^n 1/k^power_fac
        """
        return np.sum(1.0 / np.power(np.arange(1, n + 1), power_fac))

    def compute_mean_S(self, N_points, method="chi2"):
        """
        Mean of S = sum_i -log(Y_i), where Y are spacings from N_points uniforms.

        Preserves the original method switch:
        - method='chi2' uses k = N_points + 1;
        - method='gamma' uses k = N_points for the restricted statistic.
        """
        if method == "chi2":
            k = N_points + 1
        elif method == "gamma":
            k = N_points
        else:
            raise ValueError("compute_mean_S function: method must be one of: 'chi2', 'gamma'.")
        return (k - 1) * self.harmonic_numbers(k)

    def compute_variance_S(self, N_points, method="chi2"):
        """
        Variance of S = -log(prod Y_i), preserving the original formulas.
        """
        if method == "chi2":
            k = N_points + 1
            H2 = self.harmonic_numbers(k - 1, power_fac=2)
            return k**2 * H2 - (k**2 - k) * (np.pi**2) / 6.0

        elif method == "gamma":
            return (
                np.power(N_points - 1, 2) * self.harmonic_numbers(N_points, power_fac=2)
                - (N_points - 1) * (N_points - 2) * (np.pi**2) / 6
            )

        else:
            raise ValueError("compute_variance_S function: method must be one of: 'chi2', 'gamma'.")

    def exaxct_n2_cdf(self, x):
        """
        Exact survival function for N_points=2.

        Original function name is preserved for compatibility, including the
        typo in "exaxct".

            S' = -log(d_2), d_2 ~ Beta(1,2)
            SF(s) = 2 exp(-s) - exp(-2s), s >= 0
        """
        x = np.asarray(x, dtype=float)
        sf = np.zeros_like(x, dtype=float)

        mask = x >= 0
        ex = np.exp(-x[mask])
        sf[mask] = 2 * ex - ex**2

        sf = np.clip(sf, 0.0, 1.0)
        return sf

    def shifted_gamma_cdf(self, x, N_points):
        """
        Shifted Gamma approximation.

        The function name is preserved from the original script, but the value
        returned is the survival function used as a P-value.
        """
        x = np.asarray(x, dtype=float)

        mu = self.compute_mean_S(N_points, method="gamma")
        v = self.compute_variance_S(N_points, method="gamma")

        loc = (N_points - 1) * np.log(N_points - 1)
        alpha = (mu - loc) ** 2 / v
        theta = v / (mu - loc)

        return gamma.sf(x, a=alpha, loc=loc, scale=theta)

    def shifted_scaled_chi2_cdf(self, x, N_points):
        """
        Shifted-scaled chi-square approximation.

        The function name is preserved from the original script, but the value
        returned is the survival function used as a P-value.
        """
        x = np.asarray(x, dtype=float)

        df = N_points - 1
        mu = self.compute_mean_S(N_points, method="chi2")
        v = self.compute_variance_S(N_points, method="chi2")

        a_scale = np.sqrt(v / (2 * df))
        a_loc = mu - df * a_scale

        y = (x - a_loc) / a_scale

        sff = np.zeros_like(x, dtype=float)
        mask = y >= 0
        sff[mask] = chi2.sf(y[mask], df=df)

        return sff

    def S_mid_sf(self, x, N_points, method="gamma"):
        """
        Survival function for the restricted statistic:
            S' = -sum_{j=2}^n log(delta_j)

        method:
            "gamma" : shifted Gamma approximation;
            "chi2"  : shifted-scaled chi-square approximation.
        """
        if method == "gamma" and N_points < 2:
            return np.nan
        elif N_points == 2:
            return self.exaxct_n2_cdf(x)

        elif method == "chi2":
            return self.shifted_scaled_chi2_cdf(x, N_points)

        elif method == "gamma":
            return self.shifted_gamma_cdf(x, N_points)

        else:
            raise ValueError("S_mid_sf function: method must be one of: 'gamma', 'chi2'.")

    def clust_stat_from_dist_vec(self, positions_vec, pos_MT_vec, clip_ends=True):
        """
        Calculate the restricted or full spacing statistic from MR-scaled positions.

        Parameters
        ----------
        positions_vec : array-like
            Cumulative missense mutational-target coordinates.
        pos_MT_vec : array-like
            Mutational target assigned to each observed position.
        clip_ends : bool
            If True, use the restricted statistic over internal spacings only.
            If False, include endpoints 0 and 1 to form the full simplex.
        """
        positions_vec = np.asarray(positions_vec, dtype=float)
        pos_MT_vec = np.asarray(pos_MT_vec, dtype=float)

        sorted_indices = np.argsort(positions_vec)
        positions_vec = positions_vec[sorted_indices]
        pos_MT_vec = pos_MT_vec[sorted_indices]

        # Bin recurrent variants at the same coordinate.
        arr_of_arrs = [[positions_vec[0]]]
        for i in range(1, len(positions_vec)):
            if np.abs(positions_vec[i] - positions_vec[i - 1]) < 1e-10:
                arr_of_arrs[-1].append(positions_vec[i])
            else:
                arr_of_arrs.append([positions_vec[i]])

        common_ind = 0
        norm_missense_MTs = []

        # Recurrence-aware continuity correction.
        # Recurrent variants at the same position are evenly spaced within the
        # mutational target assigned to that position.
        for MT_arr in arr_of_arrs:
            L = len(MT_arr)
            for i in range(L):
                norm_missense_MT = positions_vec[common_ind] - pos_MT_vec[common_ind]
                norm_missense_MT = norm_missense_MT + (i + 1) * (
                    pos_MT_vec[common_ind] / (L + 1)
                )
                norm_missense_MTs.append(norm_missense_MT)
                common_ind += 1

        if not clip_ends:
            norm_missense_MTs = [0] + norm_missense_MTs + [1]

        diff_vec = np.diff(norm_missense_MTs)

        if np.any(diff_vec < 0):
            print(f"Warning: negative elements in diff_vec: {diff_vec}")
            print("diff_vec", diff_vec)
            print("positions_vec", positions_vec)
            print("pos_MT_vec", pos_MT_vec)
            print("arr_of_arrs", arr_of_arrs)
            print("norm_missense_MTs", norm_missense_MTs)

        stat = np.sum(-np.log(diff_vec))
        return stat

    def Sprime_sff(self, positions_vec, pos_MT_vec, method="gamma", clip_n_by=2):
        """
        Calculate the P-value/SF for the restricted clustering statistic.
        """
        n = len(positions_vec)
        if n < clip_n_by:
            return np.nan

        x = self.clust_stat_from_dist_vec(positions_vec, pos_MT_vec, clip_ends=True)
        return self.S_mid_sf(x, n, method=method)

    def cdf_permut(self, positions_vec, pos_MT_vec, permutN=1e8):
        """
        Interpolate the empirical/permutation CDF table for the clustering statistic.
        """
        if not self.D_x or not self.D_CDF:
            raise ValueError(
                "Permutation clustering requested, but no Uprod distribution file "
                "was loaded. Initialize SimplexVolumeDist(Uprod_file=...)."
            )

        n = len(positions_vec)
        if n < 2:
            return np.nan

        if n not in self.D_x:
            raise ValueError(f"No permutation clustering distribution loaded for n={n}.")

        x_arr = self.D_x[n]
        x = self.clust_stat_from_dist_vec(positions_vec, pos_MT_vec, clip_ends=True)

        ind = np.searchsorted(x_arr, x)

        if ind == 0:
            return 1.0 / permutN

        elif ind == self.L:
            return 1 - 1.0 / permutN

        cdf_arr = self.D_CDF[n]

        slope = (cdf_arr[ind] - cdf_arr[ind - 1]) / (x_arr[ind] - x_arr[ind - 1])
        intercept = cdf_arr[ind] - slope * x_arr[ind]

        return slope * x + intercept

    def calc_pval(
        self,
        positions_vec,
        pos_MT_vec,
        method="gamma",
        clip_n_by=2,
        clip_clust_pvals_by=1e-10,
    ):
        """
        Calculate the clustering P-value using gamma, chi2, or permutation method.
        """
        if method == "permut":
            pval = 1 - self.cdf_permut(positions_vec, pos_MT_vec, permutN=1e8)
        elif method in ["gamma", "chi2"]:
            pval = self.Sprime_sff(
                positions_vec,
                pos_MT_vec,
                method=method,
                clip_n_by=clip_n_by,
            )
        else:
            raise ValueError("calc_pval function: method must be one of: 'gamma', 'chi2', 'permut'.")

        if not pd.isna(pval) and pval < clip_clust_pvals_by:
            pval = clip_clust_pvals_by

        return pval


def parse_float_list(value):
    """Parse a comma-separated list of floats, returning None for missing values."""
    if pd.isna(value):
        return None
    value = str(value)
    if value == "" or value == "NA":
        return None
    return np.array(list(map(float, value.split(","))))


def calc_MisClust_P(
    merged_df,
    simplex_volume_dist,
    clip_n_by=2,
    clip_clust_pvals_by=1e-10,
    method="gamma",
    suppress_clust=False,
):
    """
    Calculate missense clustering P-values.

    Adds
    ----
    Missense_Clust_Pval
    """
    validate_columns(
        merged_df.columns,
        {"missense_MTs", "missense_MTs_by_pos"},
        "reference-annotated summary-stat DataFrame",
    )

    merged_df = merged_df.copy()

    P_vals = []
    for i in range(len(merged_df)):
        if suppress_clust:
            P_vals.append(np.nan)
            continue

        norm_mut_pos_list = merged_df["missense_MTs"].iloc[i]
        norm_position_MTs = merged_df["missense_MTs_by_pos"].iloc[i]

        positions = parse_float_list(norm_mut_pos_list)
        pos_MTs = parse_float_list(norm_position_MTs)

        if positions is None or pos_MTs is None:
            P_vals.append(np.nan)
            continue

        if len(positions) != len(pos_MTs):
            raise ValueError(
                "missense_MTs and missense_MTs_by_pos have different lengths "
                f"for row {i}: {len(positions)} vs {len(pos_MTs)}"
            )

        P_val = simplex_volume_dist.calc_pval(
            positions,
            pos_MTs,
            method=method,
            clip_n_by=clip_n_by,
            clip_clust_pvals_by=clip_clust_pvals_by,
        )

        P_vals.append(P_val)

    merged_df["Missense_Clust_Pval"] = P_vals
    return merged_df


def add_association_statistics(
    merged_df,
    syn_N=None,
    score_IDs=("AM", "REVEL", "PAI"),
    na_for0s=False,
    clust_method="gamma",
    clip_n_by=2,
    clip_clust_pvals_by=1e-10,
    suppress_clust=False,
    Uprod_file=None,
):
    """
    Add all ASCEND component P-values to a reference-annotated DataFrame.

    This performs the statistics part only:
        LoF_enrich_Pval
        AM_MisEnrich_Pval
        REVEL_MisEnrich_Pval
        PAI_MisEnrich_Pval
        Missense_Clust_Pval

    It does not combine these P-values and does not perform FDR correction.
    """
    validate_columns(merged_df.columns, {"syn_varN"}, "reference-annotated summary-stat DataFrame")

    if syn_N is None:
        syn_N = merged_df["syn_varN"].sum()

    out_df = calc_LoF_enrich_P(merged_df, syn_N=syn_N, na_for0s=na_for0s)

    for score_ID in score_IDs:
        out_df = calc_MisEnrich_P(out_df, score_ID, synN=syn_N)

    simplex_volume_dist = SimplexVolumeDist(Uprod_file=Uprod_file)

    out_df = calc_MisClust_P(
        out_df,
        simplex_volume_dist=simplex_volume_dist,
        clip_n_by=clip_n_by,
        clip_clust_pvals_by=clip_clust_pvals_by,
        method=clust_method,
        suppress_clust=suppress_clust,
    )

    return out_df


def apply_statistics_from_file(
    input_file,
    output_file,
    syn_N=None,
    na_for0s=False,
    clust_method="gamma",
    clip_n_by=2,
    clip_clust_pvals_by=1e-10,
    suppress_clust=False,
    Uprod_file=None,
):
    """
    File-level wrapper for adding ASCEND component P-values.
    """
    df = pd.read_csv(input_file, sep="\t")

    df = add_association_statistics(
        df,
        syn_N=syn_N,
        na_for0s=na_for0s,
        clust_method=clust_method,
        clip_n_by=clip_n_by,
        clip_clust_pvals_by=clip_clust_pvals_by,
        suppress_clust=suppress_clust,
        Uprod_file=Uprod_file,
    )

    df.to_csv(output_file, sep="\t", index=False, na_rep="NA")
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add ASCEND component association P-values to reference-annotated summary stats."
    )
    parser.add_argument("--i", required=True, help="Reference-annotated summary-stat TSV")
    parser.add_argument("--o", required=True, help="Output TSV with component P-values")
    parser.add_argument(
        "--syn_N",
        type=float,
        default=None,
        help="Total synonymous count for cohort scaling. Defaults to sum(syn_varN).",
    )
    parser.add_argument(
        "--r",
        action="store_true",
        help="Legacy/lenient LoF processing: convert LoF P-values >= 1 to NA.",
    )
    parser.add_argument(
        "--clust_method",
        default="gamma",
        choices=["gamma", "chi2", "permut"],
        help="Method for missense clustering P-values.",
    )
    parser.add_argument(
        "--Uprod_file",
        default=None,
        help="Required only for --clust_method permut.",
    )
    parser.add_argument(
        "--m",
        type=int,
        default=2,
        help="Minimum number of missense variants for clustering P-value calculation.",
    )
    parser.add_argument(
        "--c",
        type=float,
        default=1e-10,
        help="Lower clipping threshold for missense clustering P-values.",
    )
    parser.add_argument(
        "--suppress_clust",
        action="store_true",
        help="Suppress missense clustering P-value calculation.",
    )
    args = parser.parse_args()

    apply_statistics_from_file(
        input_file=args.i,
        output_file=args.o,
        syn_N=args.syn_N,
        na_for0s=args.r,
        clust_method=args.clust_method,
        clip_n_by=args.m,
        clip_clust_pvals_by=args.c,
        suppress_clust=args.suppress_clust,
        Uprod_file=args.Uprod_file,
    )
