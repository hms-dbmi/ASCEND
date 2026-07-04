"""
ASCEND preprocessing utilities.

This module converts observed de novo variants from a VCF-like file into per-gene
observed summary statistics used by the later ASCEND pipeline steps.
"""

from __future__ import annotations
import gzip
from dataclasses import dataclass, field
import pandas as pd

# Columns used to uniquely identify a variant in both the VCF and annotation table.
VARIANT_KEY_COLUMNS = ("#CHROM", "POS", "REF", "ALT")

# Minimal columns required in the variant annotation / mutation-target table.
# The table is expected to contain one row per possible annotated variant.
REQUIRED_ANNOTATION_COLUMNS = {
    "#CHROM",
    "POS",
    "REF",
    "ALT",
    "ensembl_gene_id",
    "Synonymous-flag",
    "LOFTEE-flag",
    "AlphaMissense_MT",
    "REVEL_MT",
    "PrimateAI-3D_MT",
}

# Output schema for per-gene observed summary statistics.
GENE_SUMMARY_COLUMNS = [
    "ENS_ID",
    "AM_y",
    "REVEL_y",
    "PAI_y",
    "Lof_varN",
    "syn_varN",
    "missense_positions",
]


def open_text(path_str: str):
    """Open a plain-text or gzip-compressed file for reading as text."""
    if path_str.endswith(".gz"):
        return gzip.open(path_str, "rt")
    return open(path_str, "r")


def make_variant_id(chrom: str, pos: str, ref: str, alt: str) -> str:
    """Create the same variant key for VCF rows and annotation-table rows."""
    return f"{chrom}_{pos}_{ref}_{alt}"


def validate_columns(columns: list[str], required: set[str], source_name: str = "") -> None:
    """Raise an error if a table is missing required columns."""
    missing = required - set(columns)
    if missing:
        label = f" in {source_name}" if source_name else ""
        raise ValueError(f"Missing required columns{label}: {sorted(missing)}")


def read_vcf_variant_counts(vcf_file: str) -> dict[str, int]:
    """
    Read a VCF/VCF.GZ and count recurrences of each observed variant.

    Returns:
        Dictionary mapping variant_id -> recurrence count, where variant_id is
        given by make_variant_id (generally CHROM_POS_REF_ALT).
    """
    variant_counts: dict[str, int] = {}

    with open_text(vcf_file) as fin:
        for line in fin:
            # Skip VCF metadata and header lines.
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed VCF line with <5 fields: {line}")

            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            variant_id = make_variant_id(chrom, pos, ref, alt)
            variant_counts[variant_id] = variant_counts.get(variant_id, 0) + 1

    #print(variant_counts)
    return variant_counts


def subset_observed_variants(
    variant_counts: dict[str, int],
    annotation_file: str,
) -> pd.DataFrame:
    """
    Subset the full annotation table to variants observed in the input VCF.

    The returned table keeps the original annotation columns and replaces
    SAMPLE_RECURRENCE with the recurrence count observed in the VCF.
    """
    rows = []

    with open_text(annotation_file) as fin:
        header = fin.readline().rstrip("\n").split("\t")
        validate_columns(header, REQUIRED_ANNOTATION_COLUMNS, str(annotation_file))
        idx = {name: i for i, name in enumerate(header)}

        for line in fin:
            fields = line.rstrip("\n").split("\t")
            variant_id = make_variant_id(
                fields[idx["#CHROM"]],
                fields[idx["POS"]],
                fields[idx["REF"]],
                fields[idx["ALT"]],
            )

            recurrence = variant_counts.get(variant_id)
            if recurrence is None:
                continue

            row = dict(zip(header, fields))
            row["SAMPLE_RECURRENCE"] = recurrence
            rows.append(row)

    # Preserve columns even when no variants match, so downstream validation
    # fails less cryptically and empty inputs remain well-formed.
    return pd.DataFrame(rows, columns=header + ["SAMPLE_RECURRENCE"])


@dataclass
class GeneObservedStats:
    """Accumulator for observed variant counts and score sums in one gene."""

    ens_id: str
    lof_count: int = 0
    syn_count: int = 0
    am_score_sum: float = 0.0
    revel_score_sum: float = 0.0
    pai_score_sum: float = 0.0
    has_am: bool = False
    has_revel: bool = False
    has_pai: bool = False
    missense_positions: list[int] = field(default_factory=list)

    def update_from_variant_row(self, row: pd.Series) -> None:
        """Update this gene's observed summary statistics using one variant row."""
        recurrence = int(row["SAMPLE_RECURRENCE"])

        # Synonymous variants are used later for cohort-level mutation-rate scaling.
        if str(row["Synonymous-flag"]) == "1":
            self.syn_count += recurrence
            return

        # LoF variants contribute only to the LoF burden count.
        if str(row["LOFTEE-flag"]) == "1":
            self.lof_count += recurrence
            return

        # Remaining coding variants with available score-based mutational targets
        # are treated as missense evidence.
        position = int(row["POS"])
        has_any_missense_score = False

        if pd.notna(row["AlphaMissense_MT"]) and row["AlphaMissense_MT"] != "NA":
            self.am_score_sum += float(row["AlphaMissense_MT"]) * recurrence
            self.has_am = True
            has_any_missense_score = True

        if pd.notna(row["REVEL_MT"]) and row["REVEL_MT"] != "NA":
            self.revel_score_sum += float(row["REVEL_MT"]) * recurrence
            self.has_revel = True
            has_any_missense_score = True

        if pd.notna(row["PrimateAI-3D_MT"]) and row["PrimateAI-3D_MT"] != "NA":
            self.pai_score_sum += float(row["PrimateAI-3D_MT"]) * recurrence
            self.has_pai = True
            has_any_missense_score = True

        # Store one position per observed recurrence. This preserves recurrence
        # information for the later missense-clustering statistic.
        if has_any_missense_score:
            self.missense_positions.extend([position] * recurrence)

    def to_record(self) -> dict:
        """Convert the accumulator to one output-table row."""
        return {
            "ENS_ID": self.ens_id,
            "AM_y": self.am_score_sum if self.has_am else pd.NA,
            "REVEL_y": self.revel_score_sum if self.has_revel else pd.NA,
            "PAI_y": self.pai_score_sum if self.has_pai else pd.NA,
            "Lof_varN": self.lof_count,
            "syn_varN": self.syn_count,
            "missense_positions": (
                ",".join(map(str, self.missense_positions))
                if self.missense_positions
                else pd.NA
            ),
        }


def aggregate_gene_observed_stats(observed_variants: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate observed annotated variants into one row per gene.

    Output columns are observed quantities only. Reference expectations and
    mutation-rate-scaled coordinates are added later in reference.py.
    """
    validate_columns(
        list(observed_variants.columns),
        REQUIRED_ANNOTATION_COLUMNS,
        "observed variant table",
    )

    gene_stats: dict[str, GeneObservedStats] = {}

    for _, row in observed_variants.iterrows():
        ens_id = row["ensembl_gene_id"]
        if ens_id not in gene_stats:
            gene_stats[ens_id] = GeneObservedStats(ens_id=ens_id)
        gene_stats[ens_id].update_from_variant_row(row)

    records = [x.to_record() for x in gene_stats.values()]
    return pd.DataFrame.from_records(records, columns=GENE_SUMMARY_COLUMNS)


def preprocess_vcf(
    vcf_file: str,
    annotation_file: str,
    keep_variant_table: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """
    Convert one VCF into observed per-gene ASCEND summary statistics.

    Args:
        vcf_file: Input VCF or VCF.GZ containing observed de novo variants.
        annotation_file: Full annotated mutation-target table.
        keep_variant_table: If True, also return the observed annotated variants.

    Returns:
        gene_stats: One row per gene with observed summary statistics.
        observed_variants: Subset annotation table, or None if not requested.
    """
    variant_counts = read_vcf_variant_counts(vcf_file)
    observed_variants = subset_observed_variants(variant_counts, annotation_file)
    #print(observed_variants)
    gene_stats = aggregate_gene_observed_stats(observed_variants)

    if keep_variant_table:
        return gene_stats, observed_variants
    return gene_stats, None


def write_tsv(df: pd.DataFrame, output_file: str) -> None:
    """Write a DataFrame as a tab-separated file using NA for missing values."""
    df.to_csv(output_file, sep="\t", index=False, na_rep="NA")


def preprocess_vcf_to_files(
    vcf_file: str,
    annotation_file: str,
    output_prefix: str,
    write_variant_table: bool = True,
) -> str:
    """
    Run VCF preprocessing and write output tables.

    Writes:
        {output_prefix}_sumstats.tsv
        {output_prefix}_variants.tsv, if write_variant_table=True

    Returns:
        Path to the observed per-gene summary-statistics file.
    """
    gene_stats, observed_variants = preprocess_vcf(
        vcf_file,
        annotation_file,
        keep_variant_table=write_variant_table,
    )

    sumstats_file = f"{output_prefix}_sumstats.tsv"
    write_tsv(gene_stats, sumstats_file)

    if write_variant_table and observed_variants is not None:
        variant_file = f"{output_prefix}_variants.tsv"
        write_tsv(observed_variants, variant_file)

    return sumstats_file
