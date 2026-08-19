#!/usr/bin/env python3
"""
Aggregate Stoner-SCF n_sweep CSV outputs into a single summary DataFrame.

Expects filenames like:
    nsweep_j0.4_-0.11111111111111112_8.0.csv
    nsweep_j<J>_<CRYSTAL_FIELD>_<NELEC>.csv

Each CSV has columns: MCA, Nelec, isYZDeg

For every file, this script computes sum(|MCA|) over all rows and
collects it into a DataFrame alongside (j, crystal_field, nelec).
"""

import re
from pathlib import Path

import pandas as pd

# Regex to pull j, crystal_field, and nelec out of the filename.
# Handles negative numbers, decimals, and long float reprs.
FNAME_PATTERN = re.compile(
    r"^nsweep_j(?P<j>-?\d+(?:\.\d+)?)_"
    r"(?P<crystal_field>-?\d+(?:\.\d+)?)_"
    r"(?P<nelec>-?\d+(?:\.\d+)?)\.csv$"
)


def parse_filename(filename: str):
    """Extract (j, crystal_field, nelec) floats from a sweep filename."""
    match = FNAME_PATTERN.match(filename)
    if not match:
        return None
    return (
        float(match.group("j")),
        float(match.group("crystal_field")),
        float(match.group("nelec")),
    )


def aggregate_sweep_dir(sweep_dir: str | Path) -> pd.DataFrame:
    """Walk sweep_dir, sum |MCA| per file, return a tidy DataFrame."""
    sweep_dir = Path(sweep_dir)
    rows = []

    for csv_path in sorted(sweep_dir.glob("*.csv")):
        parsed = parse_filename(csv_path.name)
        if parsed is None:
            print(f"Skipping (name doesn't match pattern): {csv_path.name}")
            continue
        j, crystal_field, nelec = parsed

        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"Skipping (failed to read): {csv_path.name} ({e})")
            continue

        if "MCA" not in df.columns:
            print(f"Skipping (no MCA column): {csv_path.name}")
            continue

        sum_abs_mca = df["MCA"].abs().sum()
        n_rows = len(df)

        rows.append(
            {
                "j": j,
                "crystal_field": crystal_field,
                "nelec_from_filename": nelec,
                "sum_abs_mca": sum_abs_mca,
                "n_rows": n_rows,
                "filename": csv_path.name,
            }
        )

    result = pd.DataFrame(rows).sort_values(
        ["j", "crystal_field", "nelec_from_filename"]
    ).reset_index(drop=True)

    return result


if __name__ == "__main__":
    SWEEP_DIR = "/home/cmp/Documents/Github/Stoner-SCF/out/n_sweep_csv"

    summary_df = aggregate_sweep_dir(SWEEP_DIR)

    print(summary_df.to_string(index=False))

    out_path = Path(SWEEP_DIR).parent / "mca_sweep_summary.csv"
    summary_df.to_csv(out_path, index=False)
    print(f"\nSaved summary to: {out_path}")