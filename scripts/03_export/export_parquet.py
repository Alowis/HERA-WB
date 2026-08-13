"""
Convert variable-oriented CSVs to Parquet format.

Input:
    data/SHARE/timeseries_var/*.csv  — one CSV per variable (cols = catchments)

Output:
    data/SHARE/timeseries_var_parquet/*.parquet — same structure, compressed

Usage:
    python export_parquet.py

Requirements:
    pip install pandas pyarrow
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path

# --- Paths --------------------------------------------------------------------
base_dir = Path("D:/tilloal/Documents/01_Projects/RegimeShifts")
csv_dir = base_dir / "data" / "SHARE" / "timeseries_var"
pq_dir = base_dir / "data" / "SHARE" / "timeseries_var_parquet"
pq_dir.mkdir(parents=True, exist_ok=True)

# --- Determine time vector from first file ------------------------------------
csv_files = sorted(csv_dir.glob("*.csv"))

if not csv_files:
    raise FileNotFoundError(f"No CSV files found in: {csv_dir}")

# Count rows to determine temporal resolution
with open(csv_files[0], "r") as f:
    n_rows = sum(1 for _ in f) - 1  # subtract header

n_days = (pd.Timestamp("2020-12-31") - pd.Timestamp("1951-01-01")).days + 1

if n_rows == n_days:
    print("Data is DAILY.")
    time_index = pd.date_range("1951-01-01", "2020-12-31", freq="D")
elif n_rows == n_days * 4:
    print("Data is 6-HOURLY.")
    time_index = pd.date_range("1951-01-01", "2020-12-31 18:00:00", freq="6h")
else:
    print(f"Unexpected row count: {n_rows}. Assuming daily.")
    time_index = pd.date_range("1951-01-01", periods=n_rows, freq="D")

print(f"Files: {len(csv_files)} | Timesteps: {n_rows}\n")

# --- Convert each CSV to Parquet ----------------------------------------------
total_csv_mb = 0
total_pq_mb = 0

for csv_path in csv_files:
    pq_path = pq_dir / (csv_path.stem + ".parquet")

    # Skip if already converted
    if pq_path.exists():
        print(f"  SKIP (exists): {csv_path.name}")
        continue

    print(f"  Reading: {csv_path.name} ...", end="", flush=True)
    df = pd.read_csv(csv_path)

    # Insert time as first column
    df.insert(0, "time", time_index)

    # Write Parquet with zstd compression
    print(" writing parquet...", end="", flush=True)
    df.to_parquet(pq_path, engine="pyarrow", compression="zstd", index=False)

    # Report size
    csv_size = csv_path.stat().st_size / 1e6
    pq_size = pq_path.stat().st_size / 1e6
    ratio = csv_size / pq_size
    total_csv_mb += csv_size
    total_pq_mb += pq_size

    print(f" done. ({csv_size:.0f} MB -> {pq_size:.0f} MB, {ratio:.1f}x smaller)")

    del df

print(f"\nTotal: {total_csv_mb/1000:.1f} GB -> {total_pq_mb/1000:.1f} GB "
      f"({total_csv_mb/total_pq_mb:.1f}x compression)")
print(f"Output: {pq_dir}")
print("\nTo read in Python:")
print('  df = pd.read_parquet("file.parquet")')
print('  # Or specific columns:')
print('  df = pd.read_parquet("file.parquet", columns=["time", "12345"])')
print("\nTo read in R (requires arrow package, R >= 4.3):")
print('  dt <- arrow::read_parquet("file.parquet")')
