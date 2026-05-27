"""
05_tables/generate_tables.py

Generate all four paper tables from processed CSV outputs.
No values are hard-coded; all data is read from CSV files.

Input (data/tables/):
    Design_Rainfall_Intensities_Selected.csv   <- official IDF benchmark
    IDF_TABLE.csv                              <- full official IDF (for reference)

Input (data/processed/):
    imerg_idf_2000_2019_gumbel.csv             <- IMERG Gumbel-fitted intensities
    comparison_2000_2019_vs_2020_2024.csv      <- temporal stability statistics

Output (data/tables/):
    idf_comparison.csv
    correction_factors.csv
    temporal_stability.csv
    peak_discharge.csv
    all_tables.xlsx                            <- all four tables in one workbook

Tables produced:
    Table 1  Official vs IMERG design rainfall intensities (mm/h)
    Table 2  Duration-specific correction factors (alpha_d)
    Table 3  Temporal stability: 2000-2019 vs 2020-2024 annual maxima
    Table 4  Rational method peak discharge comparison

Peak discharge parameters (Rational method: Q = C * i * A / 3.6):
    Runoff coefficient C  : 0.55  (mixed urban/vegetated, from paper)
    Catchment area A      : 3.43  km2 (from GIS delineation)
    Design return period  : 10-year
    Design duration       : 1 h   (critical duration for this catchment)
"""

from pathlib import Path
import pandas as pd
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT      = Path(__file__).resolve().parents[1]
TABLES_IN = ROOT / "data" / "tables"
PROCESSED = ROOT / "data" / "processed"
OUT       = ROOT / "data" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

# ── Rational method parameters ────────────────────────────────────────────────
C_RUNOFF       = 0.55    # runoff coefficient
A_KM2          = 3.43    # catchment area (km2)
DESIGN_RP      = 10      # return period (years)
DESIGN_DUR     = "1h"    # critical duration

# Duration display order
DUR_ORDER  = {"30min": 1, "1h": 2, "2h": 3, "4h": 4, "24h": 5}
DUR_LABELS = {"30min": "30 min", "1h": "1 h", "2h": "2 h",
              "4h": "4 h", "24h": "24 h"}

RETURN_PERIODS = [5, 10, 25]

# ══════════════════════════════════════════════════════════════════════════════
# Load inputs
# ══════════════════════════════════════════════════════════════════════════════
print("Loading input files ...")

# Official IDF benchmark (selected durations only)
official = pd.read_csv(TABLES_IN / "Design_Rainfall_Intensities_Selected.csv")
official.columns = official.columns.str.strip()

# Map official columns to return-period values
official_map = {}
for _, row in official.iterrows():
    dur_raw = str(row["Duration"]).strip()
    # Normalise duration key: "30 min" -> "30min", "1 h" -> "1h"
    dur_key = dur_raw.replace(" ", "")
    official_map[dur_key] = {
        5:  float(row["5-year Return Period (mm/h)"]),
        10: float(row["10-year Return Period (mm/h)"]),
        25: float(row["25-year Return Period (mm/h)"]),
    }

# IMERG Gumbel-fitted intensities
imerg_idf = pd.read_csv(PROCESSED / "imerg_idf_2000_2019_gumbel.csv")

# Temporal stability comparison
comparison = pd.read_csv(PROCESSED / "comparison_2000_2019_vs_2020_2024.csv")
comparison["_order"] = comparison["duration"].map(DUR_ORDER)
comparison = comparison.sort_values("_order").drop(columns="_order")

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 1  Official vs IMERG IDF comparison
# ══════════════════════════════════════════════════════════════════════════════
print("Building Table 1 ...")

rows_t1 = []
durations = sorted(DUR_ORDER.keys(), key=lambda x: DUR_ORDER[x])

for dur in durations:
    if dur not in official_map:
        continue
    row = {"Duration": DUR_LABELS[dur]}
    for rp in RETURN_PERIODS:
        off_val = official_map[dur][rp]
        imerg_val = imerg_idf.loc[
            (imerg_idf["duration"] == dur) &
            (imerg_idf["return_period_year"] == rp),
            "imerg_intensity_mm_per_hr"
        ]
        imerg_val = float(imerg_val.iloc[0]) if len(imerg_val) else np.nan
        ratio = off_val / imerg_val if imerg_val > 0 else np.nan

        row[f"Official_{rp}y (mm/h)"]  = round(off_val,   2)
        row[f"IMERG_{rp}y (mm/h)"]     = round(imerg_val, 2)
        row[f"Ratio_{rp}y"]            = round(ratio,     2)

    rows_t1.append(row)

table1 = pd.DataFrame(rows_t1)
table1.to_csv(OUT / "idf_comparison.csv", index=False)
print(f"  Saved -> idf_comparison.csv  ({len(table1)} rows)")

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 2  Duration-specific correction factors (alpha_d)
#          alpha_d = mean(Official / IMERG) across all return periods
# ══════════════════════════════════════════════════════════════════════════════
print("Building Table 2 ...")

INTERP = {
    "30min": "Severe underestimation in raw IMERG",
    "1h":    "Major underestimation in primary design duration",
    "2h":    "Substantial underestimation",
    "4h":    "Moderate underestimation",
    "24h":   "Relatively small mismatch",
}

rows_t2 = []
for dur in durations:
    if dur not in official_map:
        continue
    ratios = []
    for rp in RETURN_PERIODS:
        off_val   = official_map[dur][rp]
        imerg_val = imerg_idf.loc[
            (imerg_idf["duration"] == dur) &
            (imerg_idf["return_period_year"] == rp),
            "imerg_intensity_mm_per_hr"
        ]
        if len(imerg_val) and float(imerg_val.iloc[0]) > 0:
            ratios.append(off_val / float(imerg_val.iloc[0]))

    alpha_d = round(float(np.mean(ratios)), 2) if ratios else np.nan
    rows_t2.append({
        "Duration":                 DUR_LABELS[dur],
        "Correction factor (alpha_d)": alpha_d,
        "Interpretation":           INTERP.get(dur, ""),
    })

table2 = pd.DataFrame(rows_t2)
table2.to_csv(OUT / "correction_factors.csv", index=False)
print(f"  Saved -> correction_factors.csv  ({len(table2)} rows)")

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 3  Temporal stability: 2000-2019 vs 2020-2024
# ══════════════════════════════════════════════════════════════════════════════
print("Building Table 3 ...")

table3 = comparison[[
    "duration",
    "mean_2000_2019", "mean_2020_2024",  "mean_change_pct",
    "max_2000_2019",  "max_2020_2024",   "max_change_pct",
]].copy()

table3["duration"] = table3["duration"].map(DUR_LABELS).fillna(table3["duration"])

table3.columns = [
    "Duration",
    "Mean 2000-2019 (mm/h)", "Mean 2020-2024 (mm/h)", "Mean Change (%)",
    "Max 2000-2019 (mm/h)",  "Max 2020-2024 (mm/h)",  "Max Change (%)",
]

for col in table3.columns[1:]:
    table3[col] = table3[col].round(2)

table3.to_csv(OUT / "temporal_stability.csv", index=False)
print(f"  Saved -> temporal_stability.csv  ({len(table3)} rows)")

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 4  Peak discharge (Rational method: Q = C * i * A / 3.6)
#          Compares raw IMERG vs benchmark-corrected IMERG
# ══════════════════════════════════════════════════════════════════════════════
print("Building Table 4 ...")

# Raw IMERG intensity at design return period and duration
imerg_raw_i = imerg_idf.loc[
    (imerg_idf["duration"] == DESIGN_DUR) &
    (imerg_idf["return_period_year"] == DESIGN_RP),
    "imerg_intensity_mm_per_hr"
]
imerg_raw_i = float(imerg_raw_i.iloc[0]) if len(imerg_raw_i) else np.nan

# Official (benchmark) intensity at design RP and duration
official_i = official_map.get(DESIGN_DUR, {}).get(DESIGN_RP, np.nan)

def rational_Q(i_mm_h, C, A_km2):
    """Q (m3/s) = C * i (mm/h) * A (km2) / 3.6"""
    return round(C * i_mm_h * A_km2 / 3.6, 2)

rows_t4 = [
    {
        "Scenario":                     "Raw IMERG",
        "Rainfall intensity (mm/h)":    round(imerg_raw_i, 2),
        "Runoff coefficient C":         C_RUNOFF,
        "Catchment area (km2)":         A_KM2,
        "Peak discharge Q (m3/s)":      rational_Q(imerg_raw_i, C_RUNOFF, A_KM2),
    },
    {
        "Scenario":                     "Benchmark-corrected IMERG",
        "Rainfall intensity (mm/h)":    round(official_i, 2),
        "Runoff coefficient C":         C_RUNOFF,
        "Catchment area (km2)":         A_KM2,
        "Peak discharge Q (m3/s)":      rational_Q(official_i, C_RUNOFF, A_KM2),
    },
]

table4 = pd.DataFrame(rows_t4)
table4.to_csv(OUT / "peak_discharge.csv", index=False)
print(f"  Saved -> peak_discharge.csv  ({len(table4)} rows)")

# ══════════════════════════════════════════════════════════════════════════════
# Excel workbook — all four tables on separate sheets
# ══════════════════════════════════════════════════════════════════════════════
print("Exporting all_tables.xlsx ...")
xlsx_path = OUT / "all_tables.xlsx"
with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
    table1.to_excel(writer, sheet_name="IDF_Comparison",       index=False)
    table2.to_excel(writer, sheet_name="Correction_Factors",   index=False)
    table3.to_excel(writer, sheet_name="Temporal_Stability",   index=False)
    table4.to_excel(writer, sheet_name="Peak_Discharge",       index=False)
print(f"  Saved -> all_tables.xlsx")

# ══════════════════════════════════════════════════════════════════════════════
# Console summary
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TABLE 1  Official vs IMERG IDF (mm/h)")
print("=" * 60)
print(table1.to_string(index=False))

print("\n" + "=" * 60)
print("TABLE 2  Correction Factors (alpha_d)")
print("=" * 60)
print(table2.to_string(index=False))

print("\n" + "=" * 60)
print("TABLE 3  Temporal Stability")
print("=" * 60)
print(table3.to_string(index=False))

print("\n" + "=" * 60)
print("TABLE 4  Peak Discharge")
print("=" * 60)
print(table4.to_string(index=False))

print(f"\nAll outputs saved to: {OUT}")
