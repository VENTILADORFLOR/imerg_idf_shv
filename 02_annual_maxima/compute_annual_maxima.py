"""
02_annual_maxima/compute_annual_maxima.py

Read IMERG half-hourly series for both periods, compute rolling annual maxima,
and export all outputs needed by downstream scripts and the paper.

Input:
    data/raw/rainfall_point_imerg_2000_2019.csv
    data/raw/rainfall_point_imerg_2020_2024.csv

Output (data/processed/):
    annual_maxima_imerg_2000_2019.csv       ← fed to 03_idf_fitting/gumbel_idf_2000_2019.py
    annual_maxima_imerg_2020_2024.csv
    summary_2000_2019.csv                   ← paper Table 3 (2000-2019 side)
    summary_2020_2024.csv                   ← paper Table 3 (2020-2024 side)
    comparison_2000_2019_vs_2020_2024.csv   ← paper Table 3 (change %)
    IMERG_temporal_stability_comparison.xlsx
"""

from pathlib import Path
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT   = Path(__file__).resolve().parents[1]
RAW    = ROOT / "data" / "raw"
OUTPUT = ROOT / "data" / "processed"
OUTPUT.mkdir(parents=True, exist_ok=True)

# ── Duration definitions ──────────────────────────────────────────────────────
DURATIONS = {
    "30min": {"steps": 1,  "hours": 0.5},
    "1h":    {"steps": 2,  "hours": 1.0},
    "2h":    {"steps": 4,  "hours": 2.0},
    "4h":    {"steps": 8,  "hours": 4.0},
    "24h":   {"steps": 48, "hours": 24.0},
}
DURATION_ORDER = {"30min": 1, "1h": 2, "2h": 3, "4h": 4, "24h": 5}


# ── Core function ─────────────────────────────────────────────────────────────
def compute_annual_maxima(csv_path: Path, label: str):
    """
    Read IMERG half-hourly CSV and return:
        annual_df  : per-year annual maxima (depth & intensity) for each duration
        summary_df : mean / max / min / std intensity per duration
    """
    df = pd.read_csv(csv_path)
    df["datetime"] = pd.to_datetime(df["datetime"], errors="coerce")
    df["gpm_precipitation"] = pd.to_numeric(df["gpm_precipitation"], errors="coerce")
    df = df.dropna(subset=["datetime", "gpm_precipitation"]).copy()
    df = df.sort_values("datetime")
    df = df[~df["datetime"].duplicated(keep="first")]
    df = df.set_index("datetime").resample("30min").mean(numeric_only=True)

    # IMERG precipitation rate (mm/hr) → 30-min depth (mm)
    df["rain_mm_30min"] = df["gpm_precipitation"] * 0.5

    records = []
    for dur_name, cfg in DURATIONS.items():
        rolling_depth = (
            df["rain_mm_30min"]
            .rolling(window=cfg["steps"], min_periods=cfg["steps"])
            .sum()
        )
        temp = rolling_depth.dropna().rename("rolling_depth_mm").to_frame()
        temp["year"] = temp.index.year

        annual_max = temp.groupby("year")["rolling_depth_mm"].max().reset_index()
        annual_max["duration"]            = dur_name
        annual_max["duration_hours"]      = cfg["hours"]
        annual_max["intensity_mm_per_hr"] = annual_max["rolling_depth_mm"] / cfg["hours"]
        annual_max["period"]              = label
        records.append(annual_max)

    annual_df = pd.concat(records, ignore_index=True)
    annual_df["_order"] = annual_df["duration"].map(DURATION_ORDER)
    annual_df = annual_df.sort_values(["_order", "year"]).drop(columns="_order")

    summary_df = (
        annual_df.groupby("duration")["intensity_mm_per_hr"]
        .agg(["count", "mean", "max", "min", "std"])
        .reset_index()
    )
    summary_df["period"] = label
    summary_df["_order"] = summary_df["duration"].map(DURATION_ORDER)
    summary_df = summary_df.sort_values("_order").drop(columns="_order")

    return annual_df, summary_df


# ── Run both periods ──────────────────────────────────────────────────────────
annual_2000_2019, summary_2000_2019 = compute_annual_maxima(
    RAW / "rainfall_point_imerg_2000_2019.csv", "2000-2019"
)
annual_2020_2024, summary_2020_2024 = compute_annual_maxima(
    RAW / "rainfall_point_imerg_2020_2024.csv", "2020-2024"
)

# ── Comparison table (paper Table 3) ─────────────────────────────────────────
comparison = pd.merge(
    summary_2000_2019[["duration", "mean", "max", "min", "std"]],
    summary_2020_2024[["duration", "mean", "max", "min", "std"]],
    on="duration",
    suffixes=("_2000_2019", "_2020_2024"),
)
comparison["mean_change_pct"] = (
    (comparison["mean_2020_2024"] - comparison["mean_2000_2019"])
    / comparison["mean_2000_2019"] * 100
)
comparison["max_change_pct"] = (
    (comparison["max_2020_2024"] - comparison["max_2000_2019"])
    / comparison["max_2000_2019"] * 100
)
comparison["_order"] = comparison["duration"].map(DURATION_ORDER)
comparison = comparison.sort_values("_order").drop(columns="_order")

# ── Export CSVs ───────────────────────────────────────────────────────────────
annual_2000_2019.to_csv(OUTPUT / "annual_maxima_imerg_2000_2019.csv", index=False)
annual_2020_2024.to_csv(OUTPUT / "annual_maxima_imerg_2020_2024.csv", index=False)
summary_2000_2019.to_csv(OUTPUT / "summary_2000_2019.csv",            index=False)
summary_2020_2024.to_csv(OUTPUT / "summary_2020_2024.csv",            index=False)
comparison.to_csv(OUTPUT / "comparison_2000_2019_vs_2020_2024.csv",   index=False)

# ── Export Excel ──────────────────────────────────────────────────────────────
with pd.ExcelWriter(OUTPUT / "IMERG_temporal_stability_comparison.xlsx") as writer:
    annual_2000_2019.to_excel(writer, sheet_name="annual_max_2000_2019", index=False)
    annual_2020_2024.to_excel(writer, sheet_name="annual_max_2020_2024", index=False)
    summary_2000_2019.to_excel(writer, sheet_name="summary_2000_2019",   index=False)
    summary_2020_2024.to_excel(writer, sheet_name="summary_2020_2024",   index=False)
    comparison.to_excel(writer,        sheet_name="comparison",           index=False)

# ── Print summary ─────────────────────────────────────────────────────────────
print("\n=== 2000–2019 Summary ===")
print(summary_2000_2019.round(2).to_string(index=False))

print("\n=== 2020–2024 Summary ===")
print(summary_2020_2024.round(2).to_string(index=False))

print("\n=== Comparison (Table 3) ===")
print(comparison.round(2).to_string(index=False))

print(f"\nAll outputs saved to: {OUTPUT}")
