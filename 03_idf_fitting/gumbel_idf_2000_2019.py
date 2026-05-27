"""
03_idf_fitting/gumbel_idf_2000_2019.py

Fit a Gumbel distribution to the 2000-2019 annual maxima and derive
design intensities for return periods 5, 10, and 25 years.

Input:
    data/processed/annual_maxima_imerg_2000_2019.csv

Output (data/processed/):
    imerg_idf_2000_2019_gumbel.csv
"""

from pathlib import Path
import pandas as pd
from scipy.stats import gumbel_r

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT   = Path(__file__).resolve().parents[1]
INPUT  = ROOT / "data" / "processed" / "annual_maxima_imerg_2000_2019.csv"
OUTPUT = ROOT / "data" / "processed"
OUTPUT.mkdir(parents=True, exist_ok=True)

# ── Return periods ────────────────────────────────────────────────────────────
RETURN_PERIODS = [5, 10, 25]

# ── Load annual maxima ────────────────────────────────────────────────────────
df = pd.read_csv(INPUT)

# ── Gumbel fitting ────────────────────────────────────────────────────────────
records = []

for dur in ["30min", "1h", "2h", "4h", "24h"]:
    values = df[df["duration"] == dur]["intensity_mm_per_hr"].dropna().values

    if len(values) < 5:
        print(f"WARNING: only {len(values)} data points for {dur}, skipping.")
        continue

    loc, scale = gumbel_r.fit(values)

    for T in RETURN_PERIODS:
        p = 1 - 1 / T
        intensity = gumbel_r.ppf(p, loc=loc, scale=scale)

        records.append({
            "duration":              dur,
            "return_period_year":    T,
            "gumbel_loc":            round(loc, 4),
            "gumbel_scale":          round(scale, 4),
            "imerg_intensity_mm_per_hr": round(intensity, 4),
        })

idf_df = pd.DataFrame(records)

# ── Export ────────────────────────────────────────────────────────────────────
out_path = OUTPUT / "imerg_idf_2000_2019_gumbel.csv"
idf_df.to_csv(out_path, index=False)
print(f"Saved → {out_path}")
print(idf_df.to_string(index=False))
