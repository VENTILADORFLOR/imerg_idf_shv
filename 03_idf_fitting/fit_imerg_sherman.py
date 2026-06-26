"""
03_idf_fitting/fit_imerg_sherman.py

Fit the Sherman (1931) intensity-duration formula to the IMERG-derived
design intensities and report the parameters quoted in the manuscript
(Section 4.5 / Discussion: the IMERG temporal-structure parameters b and n).

Sherman formula:
    i = a / (t + b)^n

    where i is design rainfall intensity (mm/h), t is duration (hours),
    and a, b, n are fitted parameters. b acts as a minimum effective
    duration; n controls the rate of intensification with shorter duration.

Two reporting levels are produced:

    1. Mean level — a single (a, b, n) obtained by first averaging the
       IMERG design intensity across the 5-, 10-, and 25-year return
       periods at each duration, then fitting the Sherman formula once.
       This is the value reported as the IMERG "mean-level" parameters.

    2. Per return period — an independent (a, b, n) fitted to each of the
       5-, 10-, and 25-year IDF curves, giving the ranges of b and n
       across return periods.

The official benchmark parameters (b = 0.23, n = 0.86) are stated in the
source IDF document and are printed here only for side-by-side reference;
they are not refitted.

Input:
    data/processed/imerg_idf_2000_2019_gumbel.csv
        (from gumbel_idf_2000_2019.py; columns include
         duration, return_period_year, imerg_intensity_mm_per_hr)

Output (data/processed/):
    imerg_sherman_parameters.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT   = Path(__file__).resolve().parents[1]
INPUT  = ROOT / "data" / "processed" / "imerg_idf_2000_2019_gumbel.csv"
OUTPUT = ROOT / "data" / "processed"
OUTPUT.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUTPUT / "imerg_sherman_parameters.csv"

# ── Configuration ─────────────────────────────────────────────────────────────
DURATION_HOURS = {
    "30min": 0.5,
    "1h":    1.0,
    "2h":    2.0,
    "4h":    4.0,
    "24h":   24.0,
}
RETURN_PERIODS    = [5, 10, 25]
INTENSITY_COL     = "imerg_intensity_mm_per_hr"
SHERMAN_P0        = [100.0, 1.0, 0.7]   # initial guess for [a, b, n]
SHERMAN_MAXFEV    = 20000

# Official benchmark Sherman parameters (stated in the source IDF document)
OFFICIAL_B = 0.23
OFFICIAL_N = 0.86

# ── Sherman model ─────────────────────────────────────────────────────────────
def sherman(t, a, b, n):
    """Sherman intensity-duration formula: i = a / (t + b)^n."""
    return a / np.power(t + b, n)

def fit_sherman(t, i):
    """Fit Sherman parameters (a, b, n) to duration/intensity arrays."""
    popt, _ = curve_fit(sherman, t, i, p0=SHERMAN_P0, maxfev=SHERMAN_MAXFEV)
    pred = sherman(np.asarray(t), *popt)
    ss_res = np.sum((np.asarray(i) - pred) ** 2)
    ss_tot = np.sum((np.asarray(i) - np.mean(i)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return popt, r2

# ── Load IMERG IDF intensities ────────────────────────────────────────────────
df = pd.read_csv(INPUT)
df = df[df["duration"].isin(DURATION_HOURS)].copy()
df["t_h"] = df["duration"].map(DURATION_HOURS)

records = []

# ── 1. Mean level: average intensity across return periods, then fit once ─────
mean_i = (
    df[df["return_period_year"].isin(RETURN_PERIODS)]
    .groupby("t_h")[INTENSITY_COL]
    .mean()
    .sort_index()
)
(a_m, b_m, n_m), r2_m = fit_sherman(mean_i.index.values, mean_i.values)
records.append({
    "level":      "mean",
    "a":          round(a_m, 3),
    "b":          round(b_m, 3),
    "n":          round(n_m, 3),
    "r2":         round(r2_m, 5),
})

# ── 2. Per return period ──────────────────────────────────────────────────────
for T in RETURN_PERIODS:
    sub = df[df["return_period_year"] == T].sort_values("t_h")
    (a_t, b_t, n_t), r2_t = fit_sherman(sub["t_h"].values, sub[INTENSITY_COL].values)
    records.append({
        "level":  f"{T}yr",
        "a":      round(a_t, 3),
        "b":      round(b_t, 3),
        "n":      round(n_t, 3),
        "r2":     round(r2_t, 5),
    })

out_df = pd.DataFrame(records)
out_df.to_csv(OUT_CSV, index=False)

# ── Report ────────────────────────────────────────────────────────────────────
rp_rows = out_df[out_df["level"].str.endswith("yr")]
b_lo, b_hi = rp_rows["b"].min(), rp_rows["b"].max()
n_lo, n_hi = rp_rows["n"].min(), rp_rows["n"].max()

print("IMERG Sherman parameters  (i = a / (t + b)^n)")
print(out_df.to_string(index=False))
print()
print(f"Mean level:           b = {b_m:.2f}, n = {n_m:.2f}")
print(f"Across return periods: b = {b_lo:.2f}-{b_hi:.2f}, n = {n_lo:.2f}-{n_hi:.2f}")
print(f"Official benchmark:    b = {OFFICIAL_B:.2f}, n = {OFFICIAL_N:.2f} (stated in source)")
print(f"\nSaved -> {OUT_CSV}")
