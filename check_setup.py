#!/usr/bin/env python3

from pathlib import Path
import importlib
import sys

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"

REQUIRED_FILES = [
    "dominant_genes_ENS.txt",
    "gene_wFDR_weights_Sfacs50_quant.tsv",
    "Uprod_dists.txt.gz",
    "ALLVARS_MR_dist_by_pos5.txt.gz",
    "ENS_ID2Gene_ID.txt.gz",
    "BY_GENE_MR_5.txt.gz",
    "muttargs5_noOL_RQC_format.txt.gz",
]

REQUIRED_PACKAGES = [
    "numpy",
    "pandas",
    "scipy",
]


def human_size(n_bytes):
    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(n_bytes)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.1f} {unit}"
        size /= 1024


def main():
    failed = False

    print("ASCEND setup check")
    print()

    print("Python:")
    print(f"  {sys.version.split()[0]}")
    if sys.version_info < (3, 10):
        print("  ERROR: Python >=3.10 is recommended.")
        failed = True
    else:
        print("  OK")
    print()

    print("Python packages:")
    for pkg in REQUIRED_PACKAGES:
        try:
            mod = importlib.import_module(pkg)
            version = getattr(mod, "__version__", "unknown")
            print(f"  {pkg}: {version} OK")
        except ImportError:
            print(f"  {pkg}: MISSING")
            failed = True
    print()

    print(f"Reference files in {DATA_DIR}:")
    for filename in REQUIRED_FILES:
        path = DATA_DIR / filename
        if path.exists():
            print(f"  {filename}: OK ({human_size(path.stat().st_size)})")
        else:
            print(f"  {filename}: MISSING")
            failed = True
    print()

    if failed:
        print("Setup check failed. Install missing packages and/or download the ASCEND reference files.")
        sys.exit(1)

    print("Setup looks complete.")
    print()
    print("Try the toy example:")
    print("  cd example")
    print("  bash lazy-test.sh")


if __name__ == "__main__":
    main()
