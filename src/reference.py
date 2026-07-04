"""
reference.py

Utilities for adding ASCEND reference annotations to observed per-gene
summary statistics.

This module starts from the output of preprocess.py, which contains only
observed cohort data, and adds the reference information needed downstream:

1. all processed genes including genes with no observed
   variants in the cohort;
2. gene-level mutation rate expectations for LoF and missense classes;
3. common gene names;
4. mutation rate-scaled coordinates for observed missense positions.
"""

from __future__ import annotations
import gzip
import numpy as np
import pandas as pd


# Output of preprocess.py
OBSERVED_COLUMNS = {
    "ENS_ID",
    "AM_y",
    "REVEL_y",
    "PAI_y",
    "Lof_varN",
    "syn_varN",
    "missense_positions",
}

# Minimal reference table expected after column-name standardization.
REFERENCE_COLUMNS = {
    "ENS_ID",
    "AM_lambda",
    "REVEL_lambda",
    "PAI_lambda",
    "Lof_lambda",
    "syn_lambda",
    "missense_lambda",
}

# Column-name mapping from the current mutation-rate reference files to the
# standardized names used by ASCEND downstream.
REFERENCE_RENAME_MAP = {
    "Gene_ID": "ENS_ID",
    "AlphaMissense": "AM_lambda",
    "REVEL": "REVEL_lambda",
    "PrimateAI-3D": "PAI_lambda",
    "LOFTEE": "Lof_lambda",
    "SYN": "syn_lambda",
    "MissenseTotal": "missense_lambda",
}

# Legacy output names from previous versions
OBSERVED_RENAME_MAP = {
    "#ENS_ID": "ENS_ID",
    "AM": "AM_y",
    "REVEL": "REVEL_y",
    "PAI": "PAI_y",
}

OBSERVED_NUMERIC_COLUMNS = ["AM_y", "REVEL_y", "PAI_y", "Lof_varN", "syn_varN"]
OBSERVED_COUNT_COLUMNS = ["Lof_varN", "syn_varN"]
OBSERVED_SCORE_COLUMNS = ["AM_y", "REVEL_y", "PAI_y"]
LAMBDA_COLUMNS = [
    "AM_lambda",
    "REVEL_lambda",
    "PAI_lambda",
    "Lof_lambda",
    "syn_lambda",
    "missense_lambda",
]

OUTPUT_COLUMNS = [
    "ENS_ID",
    "Gene_Name",
    "AM_lambda",
    "REVEL_lambda",
    "PAI_lambda",
    "Lof_lambda",
    "syn_lambda",
    "missense_lambda",
    "AM_y",
    "REVEL_y",
    "PAI_y",
    "Lof_varN",
    "syn_varN",
    "missense_positions",
    "missense_MTs",
    "missense_MTs_by_pos",
]


def open_text(path_str: str):
    """Open a plain-text or gzipped text file."""
    if path_str.endswith(".gz"):
        return gzip.open(path_str, "rt")
    return open(path_str, "r")


def read_tsv(path_str: str, **kwargs) -> pd.DataFrame:
    """Read a tab-separated file, transparently handling .gz files."""
    with open_text(path_str) as fin:
        return pd.read_table(fin, sep="\t", **kwargs)


def write_tsv(df: pd.DataFrame, output_file: str) -> None:
    """Write a DataFrame as a tab-separated file using NA for missing values."""
    df.to_csv(output_file, sep="\t", index=False, na_rep="NA")


def validate_columns(columns: list[str], required: set[str], source_name: str = "") -> None:
    """Raise a readable error if required columns are absent."""
    missing = required - set(columns)
    if missing:
        label = f" in {source_name}" if source_name else ""
        raise ValueError(f"Missing required columns{label}: {sorted(missing)}")


def check_unique_ids(df: pd.DataFrame, id_column: str, source_name: str) -> None:
    """Make sure one gene ID appears at most once in a table."""
    duplicated = df[id_column][df[id_column].duplicated()].unique()
    if len(duplicated) > 0:
        example = ", ".join(map(str, duplicated[:5]))
        raise ValueError(f"Duplicated {id_column} values in {source_name}: {example}")


def standardize_observed_sumstats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize observed per-gene summary stats from preprocess.py.

    The preferred input columns are:
        ENS_ID, AM_y, REVEL_y, PAI_y, Lof_varN, syn_varN, missense_positions

    The old RaMeDiES names (#ENS_ID, AM, REVEL, PAI) are also accepted.
    Missing observed score/count values are converted to zero, because a gene
    with no observed variants should later contribute neutral evidence rather
    than a missing test.
    """
    df = df.copy()
    df = df.rename(columns=OBSERVED_RENAME_MAP)

    validate_columns(list(df.columns), {"ENS_ID"}, "observed summary-stat table")
    check_unique_ids(df, "ENS_ID", "observed summary-stat table")

    # Keep the module tolerant to minimal input files by creating absent observed
    # columns. This is useful for toy examples and legacy files.
    for col in OBSERVED_SCORE_COLUMNS + OBSERVED_COUNT_COLUMNS:
        if col not in df.columns:
            df[col] = 0

    if "missense_positions" not in df.columns:
        df["missense_positions"] = pd.NA

    for col in OBSERVED_NUMERIC_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in OBSERVED_COUNT_COLUMNS:
        df[col] = df[col].astype(int)

    return df[["ENS_ID"] + OBSERVED_NUMERIC_COLUMNS + ["missense_positions"]]


def read_observed_sumstats(sumstats_file: str) -> pd.DataFrame:
    """Read and standardize observed per-gene summary stats."""
    return standardize_observed_sumstats(read_tsv(sumstats_file))


def standardize_by_gene_reference(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize the by-gene mutation-rate reference table.

    The expected raw reference columns are currently:
        Gene_ID, AlphaMissense, REVEL, PrimateAI-3D, LOFTEE, SYN, MissenseTotal

    They are renamed to standardized ASCEND names ending in _lambda. These
    lambdas are assumed to be on the per-observed-synonymous-variant scale;
    downstream statistics multiply them by the cohort's observed synonymous
    count.
    """
    df = df.copy()
    df = df.rename(columns=REFERENCE_RENAME_MAP)

    validate_columns(list(df.columns), REFERENCE_COLUMNS, "by-gene mutation-rate reference")
    check_unique_ids(df, "ENS_ID", "by-gene mutation-rate reference")

    for col in LAMBDA_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def read_by_gene_reference(by_gene_mr_file: str) -> pd.DataFrame:
    """Read and standardize the by-gene mutation-rate reference table."""
    return standardize_by_gene_reference(read_tsv(by_gene_mr_file))


def merge_observed_with_reference(
    observed_df: pd.DataFrame,
    reference_df: pd.DataFrame,
    fail_on_missing_reference_genes: bool = True,
) -> pd.DataFrame:
    """
    Merge observed cohort summary stats onto the complete gene reference.

    The merge is a left join from the reference table. This adds genes with no
    observed variants in the cohort, while preserving all genes that can be
    tested by ASCEND.
    """
    observed_df = standardize_observed_sumstats(observed_df)
    reference_df = standardize_by_gene_reference(reference_df)

    if fail_on_missing_reference_genes:
        missing_genes = set(observed_df["ENS_ID"]) - set(reference_df["ENS_ID"])
        if missing_genes:
            example = ", ".join(sorted(missing_genes)[:5])
            raise ValueError(
                "Observed genes absent from the by-gene mutation-rate reference: "
                f"{example}"
            )

    merged_df = pd.merge(reference_df, observed_df, on="ENS_ID", how="left", validate="one_to_one")

    # After the reference-left merge, genes with no observed variants have NA in
    # observed columns. Convert those to zeros so downstream P-value functions can
    # return neutral P values instead of treating the tests as missing.
    for col in OBSERVED_SCORE_COLUMNS + OBSERVED_COUNT_COLUMNS:
        merged_df[col] = pd.to_numeric(merged_df[col], errors="coerce")

    if "missense_positions" not in merged_df.columns:
        merged_df["missense_positions"] = pd.NA

    return merged_df


def read_gene_names(gene_name_file: str) -> pd.DataFrame:
    """
    Read an ENS_ID-to-gene-symbol mapping file.

    The current ASCEND/RaMeDiES mapping file is headerless with two columns:
        ENS_ID    Gene_Name

    A header row with these same names is tolerated and dropped.
    """
    df = read_tsv(gene_name_file, header=None)
    if df.shape[1] < 2:
        raise ValueError("Gene-name file must contain at least two columns: ENS_ID and Gene_Name")

    df = df.iloc[:, :2].copy()
    df.columns = ["ENS_ID", "Gene_Name"]

    # Tolerate accidental header rows in otherwise headerless files.
    header_mask = (df["ENS_ID"] == "ENS_ID") & (df["Gene_Name"] == "Gene_Name")
    df = df.loc[~header_mask].copy()

    df = df.drop_duplicates(subset=["ENS_ID"], keep="first")
    return df


def add_gene_names(merged_df: pd.DataFrame, gene_name_file: str) -> pd.DataFrame:
    """Add common gene names to a reference-merged ASCEND table."""
    gene_name_df = read_gene_names(gene_name_file)

    df = merged_df.copy()
    if "Gene_Name" in df.columns:
        df = df.drop(columns=["Gene_Name"])

    df = pd.merge(df, gene_name_df, on="ENS_ID", how="left", validate="one_to_one")
    return move_column_after(df, "Gene_Name", "ENS_ID")


def move_column_after(df: pd.DataFrame, column: str, after: str) -> pd.DataFrame:
    """Return a copy of df with one column moved immediately after another."""
    if column not in df.columns or after not in df.columns:
        return df

    cols = [c for c in df.columns if c != column]
    insert_at = cols.index(after) + 1
    cols.insert(insert_at, column)
    return df[cols]


def is_missing_value(value) -> bool:
    """Treat pandas NA, empty strings, and literal NA/nan strings as missing."""
    if pd.isna(value):
        return True
    if isinstance(value, str) and value.strip().lower() in {"", "na", "nan", "none"}:
        return True
    return False


def parse_int_list(value) -> list[int]:
    """Parse comma-separated integer-like values from a TSV cell."""
    if is_missing_value(value):
        return []

    if isinstance(value, (int, np.integer)):
        return [int(value)]

    if isinstance(value, (float, np.floating)):
        return [int(value)]

    return [int(float(x)) for x in str(value).split(",") if x.strip()]


def format_float_list(values: list[float]) -> str:
    """Format floats compactly for comma-separated TSV cells."""
    if len(values) == 0:
        return pd.NA
    return ",".join(f"{x:.12g}" for x in values)


def collect_requested_positions(merged_df: pd.DataFrame) -> dict[str, set[int]]:
    """
    Collect the observed missense positions for which MR coordinates are needed.

    Returns:
        {ENS_ID: {position_1, position_2, ...}}
    """
    requested: dict[str, set[int]] = {}

    for _, row in merged_df.iterrows():
        positions = parse_int_list(row.get("missense_positions", pd.NA))
        if not positions:
            continue

        gene = row["ENS_ID"]
        requested.setdefault(gene, set()).update(positions)

    return requested


def infer_mr_distance_columns(header: list[str]) -> tuple[int, int, int]:
    """
    Infer gene, position, and cumulative-MR columns in the MR-distance file.

    The current file is expected to have these as the first three columns. This
    helper also tolerates common explicit column names.
    """
    name_to_idx = {name: idx for idx, name in enumerate(header)}

    gene_candidates = ["ENS_ID", "Gene_ID", "ensembl_gene_id", "gene"]
    pos_candidates = ["POS", "position", "pos"]
    mt_candidates = ["MT", "mutational_target", "missense_MT", "cumulative_MT"]

    gene_idx = next((name_to_idx[x] for x in gene_candidates if x in name_to_idx), 0)
    pos_idx = next((name_to_idx[x] for x in pos_candidates if x in name_to_idx), 1)
    mt_idx = next((name_to_idx[x] for x in mt_candidates if x in name_to_idx), 2)

    return gene_idx, pos_idx, mt_idx


def read_mr_coordinate_lookup(
    mr_dist_file: str,
    requested_positions: dict[str, set[int]],
) -> dict[str, dict[int, tuple[float, float]]]:
    """
    Read mutation-rate-scaled coordinates for requested missense positions.

    The MR-distance file is expected to be sorted by gene and genomic position.
    For each position, it should contain a cumulative mutation-rate coordinate
    within the gene. This function returns both:

        1. the cumulative coordinate at the right edge of the position's bin;
        2. the mutation-rate width of the position's bin.

    The clustering code then applies a recurrence-aware continuity correction
    within each position-specific bin.
    """
    lookup: dict[str, dict[int, tuple[float, float]]] = {
        gene: {} for gene in requested_positions
    }

    if not requested_positions:
        return lookup

    with open_text(mr_dist_file) as fin:
        header = fin.readline().rstrip("\n").split("\t")
        gene_idx, pos_idx, mt_idx = infer_mr_distance_columns(header)

        previous_gene = None
        previous_mt = 0.0

        for line_number, line in enumerate(fin, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(gene_idx, pos_idx, mt_idx):
                raise ValueError(f"Malformed MR-distance line {line_number}: {line[:100]}")

            gene = fields[gene_idx]
            pos = int(float(fields[pos_idx]))
            cumulative_mt = float(fields[mt_idx])

            if gene != previous_gene:
                previous_gene = gene
                previous_mt = 0.0

            position_mt = cumulative_mt - previous_mt
            if position_mt < -1e-12:
                raise ValueError(
                    "MR-distance file appears unsorted or cumulative coordinates decreased: "
                    f"gene={gene}, pos={pos}, line={line_number}"
                )

            # Clip tiny negative values caused by floating-point noise.
            position_mt = max(position_mt, 0.0)

            if gene in requested_positions and pos in requested_positions[gene]:
                lookup[gene][pos] = (cumulative_mt, position_mt)

            previous_mt = cumulative_mt

    return lookup


def add_missense_mr_coordinates(
    merged_df: pd.DataFrame,
    mr_dist_file: str,
    missing_position: str = "raise",
) -> pd.DataFrame:
    """
    Add mutation-rate-scaled coordinates for observed missense positions.

    Adds two comma-separated columns aligned with missense_positions:

        missense_MTs:
            cumulative mutation-rate coordinate for the observed position;
        missense_MTs_by_pos:
            mutation-rate width of the observed position's bin.

    Args:
        merged_df:
            Output of merge_observed_with_reference(), optionally already with
            gene names added.
        mr_dist_file:
            Position-level mutation-rate coordinate file.
        missing_position:
            "raise" to fail if an observed position is absent from the MR file;
            "ignore" to write NA coordinates for missing positions.
    """
    if missing_position not in {"raise", "ignore"}:
        raise ValueError("missing_position must be either 'raise' or 'ignore'")

    df = merged_df.copy()
    requested_positions = collect_requested_positions(df)
    coordinate_lookup = read_mr_coordinate_lookup(mr_dist_file, requested_positions)

    mts_per_gene = []
    mts_per_pos = []

    for _, row in df.iterrows():
        gene = row["ENS_ID"]
        positions = parse_int_list(row.get("missense_positions", pd.NA))

        if not positions:
            mts_per_gene.append(pd.NA)
            mts_per_pos.append(pd.NA)
            continue

        cumulative_values = []
        position_widths = []

        for pos in positions:
            if gene not in coordinate_lookup or pos not in coordinate_lookup[gene]:
                if missing_position == "raise":
                    raise KeyError(
                        "Observed missense position absent from MR-distance file: "
                        f"ENS_ID={gene}, POS={pos}"
                    )
                cumulative_values.append(np.nan)
                position_widths.append(np.nan)
                continue

            cumulative_mt, position_mt = coordinate_lookup[gene][pos]
            cumulative_values.append(cumulative_mt)
            position_widths.append(position_mt)

        mts_per_gene.append(format_float_list(cumulative_values))
        mts_per_pos.append(format_float_list(position_widths))

    df["missense_MTs"] = mts_per_gene
    df["missense_MTs_by_pos"] = mts_per_pos
    return df


def reorder_reference_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Put standard ASCEND columns first and preserve any extra columns after."""
    standard = [col for col in OUTPUT_COLUMNS if col in df.columns]
    extra = [col for col in df.columns if col not in standard]
    return df[standard + extra]


def add_reference_annotations(
    observed_sumstats: pd.DataFrame,
    by_gene_mr_file: str,
    gene_name_file: str,
    mr_dist_file: str,
    missing_position: str = "raise",
) -> pd.DataFrame:
    """
    Add all reference annotations needed for ASCEND analysis.

    Args:
        observed_sumstats:
            DataFrame from preprocess.py, with one row per observed gene.
        by_gene_mr_file:
            Gene-level mutation-rate reference file.
        gene_name_file:
            ENS_ID-to-gene-symbol mapping file.
        mr_dist_file:
            Position-level mutation-rate coordinate file.
        missing_position:
            Passed to add_missense_mr_coordinates().

    Returns:
        A complete gene-level table ready for statistics.py.
    """
    reference_df = read_by_gene_reference(by_gene_mr_file)
    merged_df = merge_observed_with_reference(observed_sumstats, reference_df)
    merged_df = add_gene_names(merged_df, gene_name_file)
    merged_df = add_missense_mr_coordinates(
        merged_df,
        mr_dist_file,
        missing_position=missing_position,
    )
    return reorder_reference_columns(merged_df)


def update_sumstats_file_with_reference(
    observed_sumstats_file: str,
    by_gene_mr_file: str,
    gene_name_file: str,
    mr_dist_file: str,
    output_file: str,
    missing_position: str = "raise",
) -> str:
    """
    File-based wrapper around add_reference_annotations().

    This is the replacement for the old sum_stats_update.py workflow, but with
    all reference paths supplied explicitly instead of hard-coded module globals.
    """
    observed_sumstats = read_observed_sumstats(observed_sumstats_file)
    annotated_df = add_reference_annotations(
        observed_sumstats=observed_sumstats,
        by_gene_mr_file=by_gene_mr_file,
        gene_name_file=gene_name_file,
        mr_dist_file=mr_dist_file,
        missing_position=missing_position,
    )
    write_tsv(annotated_df, output_file)
    return output_file
