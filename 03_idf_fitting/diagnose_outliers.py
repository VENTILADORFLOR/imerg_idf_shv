"""
diagnose_outliers.py

Local diagnostic script to identify extreme annual maxima events
in IMERG (2000-2024) and CHIRPS (1981-2025) at the Sihanoukville study point.

Run from the 03_idf_fitting/ directory:
    python diagnose_outliers.py

Inputs:
    data/raw/rainfall_point_imerg_2000_2019.csv
    data/raw/rainfall_point_imerg_2020_2024.csv
    data/raw/rainfall_point_chirps_1981_2025.csv

Outputs (data/processed/):
    outliers_imerg_annual_maxima.csv
    outliers_chirps_annual_maxima.csv
    outliers_crosscheck.csv
    outliers_log.txt
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from io import StringIO

# ── Paths ─────────────────────────────────────────────────────────────────────
# Script lives in 03_idf_fitting/ — go up one level to project root
ROOT      = Path(__file__).resolve().parent.parent
RAW       = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

IMERG_A   = RAW / "rainfall_point_imerg_2000_2019.csv"
IMERG_B   = RAW / "rainfall_point_imerg_2020_2024.csv"
CHIRPS    = RAW / "rainfall_point_chirps_1981_2025.csv"

OUT_IMERG      = PROCESSED / "outliers_imerg_annual_maxima.csv"
OUT_CHIRPS     = PROCESSED / "outliers_chirps_annual_maxima.csv"
OUT_CROSSCHECK = PROCESSED / "outliers_crosscheck.csv"
OUT_LOG        = PROCESSED / "outliers_log.txt"

# ── Tee: write to both console and log buffer ─────────────────────────────────
log_buffer = StringIO()

def log(msg=""):
    print(msg)
    log_buffer.write(msg + "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 1. IMERG — find date of each year's 24h maximum
# ══════════════════════════════════════════════════════════════════════════════
log("=" * 60)
log("IMERG 24h annual maxima — top 10 events")
log("=" * 60)

df_a = pd.read_csv(IMERG_A)
df_b = pd.read_csv(IMERG_B)
df_imerg = pd.concat([df_a, df_b], ignore_index=True)
df_imerg["datetime"] = pd.to_datetime(df_imerg["datetime"])
df_imerg = df_imerg.dropna(subset=["gpm_precipitation"])
df_imerg = df_imerg.drop_duplicates(subset=["datetime"])
df_imerg = df_imerg.sort_values("datetime").set_index("datetime")
df_imerg = df_imerg.resample("30min").mean(numeric_only=True)
df_imerg["depth_mm"] = df_imerg["gpm_precipitation"] * 0.5
df_imerg["rolling_24h_mm"] = df_imerg["depth_mm"].rolling(48).sum()
df_imerg["year"] = df_imerg.index.year

imerg_records = []
for year, grp in df_imerg.groupby("year"):
    idx_max = grp["rolling_24h_mm"].idxmax()
    val     = grp["rolling_24h_mm"].max()
    if pd.notna(val):
        imerg_records.append({
            "year":         year,
            "date_of_max":  idx_max.date(),
            "24h_depth_mm": round(val, 1),
            "rank":         None,
        })

imerg_am = pd.DataFrame(imerg_records).sort_values("24h_depth_mm", ascending=False).reset_index(drop=True)
imerg_am["rank"] = imerg_am.index + 1

log(imerg_am.head(10)[["rank","year","date_of_max","24h_depth_mm"]].to_string(index=False))
log()
log("All IMERG annual maxima (chronological):")
log(imerg_am.sort_values("year")[["year","date_of_max","24h_depth_mm","rank"]].to_string(index=False))

imerg_am.sort_values("year").to_csv(OUT_IMERG, index=False)
log(f"\nSaved → {OUT_IMERG}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. CHIRPS — find date of each year's daily maximum
# ══════════════════════════════════════════════════════════════════════════════
log()
log("=" * 60)
log("CHIRPS daily annual maxima — top 10 events")
log("=" * 60)

df_chirps = pd.read_csv(CHIRPS)
df_chirps["date"] = pd.to_datetime(df_chirps["date"])
df_chirps = df_chirps.dropna(subset=["chirps_precipitation"])
df_chirps = df_chirps[df_chirps["chirps_precipitation"] >= 0]
df_chirps["year"] = df_chirps["date"].dt.year

chirps_records = []
for year, grp in df_chirps.groupby("year"):
    idx_max  = grp["chirps_precipitation"].idxmax()
    val      = grp["chirps_precipitation"].max()
    date_max = grp.loc[idx_max, "date"].date()
    chirps_records.append({
        "year":         year,
        "date_of_max":  date_max,
        "24h_depth_mm": round(val, 1),
        "rank":         None,
    })

chirps_am = pd.DataFrame(chirps_records).sort_values("24h_depth_mm", ascending=False).reset_index(drop=True)
chirps_am["rank"] = chirps_am.index + 1

log(chirps_am.head(10)[["rank","year","date_of_max","24h_depth_mm"]].to_string(index=False))
log()
log("All CHIRPS annual maxima (chronological):")
log(chirps_am.sort_values("year")[["year","date_of_max","24h_depth_mm","rank"]].to_string(index=False))

chirps_am.sort_values("year").to_csv(OUT_CHIRPS, index=False)
log(f"\nSaved → {OUT_CHIRPS}")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Cross-check: top 10 CHIRPS events vs IMERG same year
# ══════════════════════════════════════════════════════════════════════════════
log()
log("=" * 60)
log("Cross-check: top 10 CHIRPS events vs IMERG same year")
log("=" * 60)

imerg_by_year  = imerg_am.set_index("year")
crosscheck_rows = []

log(f"{'Year':>6} {'CHIRPS date':>14} {'CHIRPS mm':>10} {'IMERG date':>14} {'IMERG mm':>10} {'Notes'}")
log("-" * 75)

for _, row in chirps_am.head(10).iterrows():
    yr     = row["year"]
    c_date = str(row["date_of_max"])
    c_val  = row["24h_depth_mm"]

    if yr in imerg_by_year.index:
        i_date = str(imerg_by_year.loc[yr, "date_of_max"])
        i_val  = float(imerg_by_year.loc[yr, "24h_depth_mm"])
        same_date = "same date" if c_date == i_date else "different date"
    else:
        i_date, i_val, same_date = "N/A", float("nan"), "no IMERG coverage"

    log(f"{yr:>6} {c_date:>14} {c_val:>10} {i_date:>14} {str(i_val):>10}  {same_date}")
    crosscheck_rows.append({
        "year":            yr,
        "chirps_date":     c_date,
        "chirps_mm":       c_val,
        "imerg_date":      i_date,
        "imerg_mm":        i_val,
        "notes":           same_date,
    })

log()
log("Note: IMERG starts 2000-06, so pre-2000 CHIRPS events have no IMERG match.")

crosscheck_df = pd.DataFrame(crosscheck_rows)
crosscheck_df.to_csv(OUT_CROSSCHECK, index=False)
log(f"\nSaved → {OUT_CROSSCHECK}")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Save log
# ══════════════════════════════════════════════════════════════════════════════
with open(OUT_LOG, "w", encoding="utf-8") as f:
    f.write(log_buffer.getvalue())
print(f"Log saved → {OUT_LOG}")

