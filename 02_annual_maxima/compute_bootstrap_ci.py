"""
02_idf_analysis/compute_bootstrap_ci.py

Compute bootstrap 95% confidence intervals for the duration-specific
IMERG correction factors (alpha_d) reported in Table 3 of the manuscript.

Method:
    For each duration, resample the annual maxima series (2000-2019) with
    replacement N_BOOT times. For each bootstrap sample, fit a Gumbel
    distribution using L-moments and extract design intensities at 5-, 10-,
    and 25-year return periods. The correction factor for each bootstrap
    replicate is the mean ratio across the three return periods. The 2.5th
    and 97.5th percentiles of the bootstrap distribution define the 95% CI.

Inputs:
    data/raw/rainfall_point_imerg_2000_2019.csv

Outputs:
    data/processed/bootstrap_ci.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT      = Path(__file__).resolve().parents[1]
CSV_PATH  = ROOT / "data" / "raw"       / "rainfall_point_imerg_2000_2019.csv"
OUT_PATH  = ROOT / "data" / "processed" / "bootstrap_ci.csv"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
DURATIONS_MIN   = [30, 60, 120, 240, 1440]
DURATIONS_LABEL = ["30min", "1h", "2h", "4h", "24h"]
RETURN_PERIODS  = [5, 10, 25]
N_BOOT          = 2000
SEED            = 42

OFFICIAL = {
    30:   {5: 144.2, 10: 172.0, 25: 207.1},
    60:   {5:  92.0, 10: 109.8, 25: 132.2},
    120:  {5:  55.2, 10:  65.8, 25:  79.3},
    240:  {5:  31.8, 10:  38.0, 25:  45.7},
    1440: {5:   7.1, 10:   8.5, 25:  10.2},
}

# ── Gumbel helpers ────────────────────────────────────────────────────────────
def lmoments_gumbel(x):
    x_sorted = np.sort(x)
    n = len(x_sorted)
    i = np.arange(1, n + 1)
    b0 = np.mean(x_sorted)
    b1 = np.sum((i - 1) / (n - 1) * x_sorted) / n
    lam2  = 2 * b1 - b0
    scale = lam2 / np.log(2)
    loc   = b0 - 0.5772 * scale
    return loc, scale

def gumbel_quantile(T, loc, scale):
    p = 1 - 1 / T
    return loc - scale * np.log(-np.log(p))

# ── Load data ─────────────────────────────────────────────────────────────────
print(f"Loading: {CSV_PATH}")
df = pd.read_csv(CSV_PATH, parse_dates=["datetime"])
df = df.rename(columns={"gpm_precipitation": "precip_mmh"}).dropna(subset=["precip_mmh"])
df["depth_mm"] = df["precip_mmh"] * 0.5
df = df[df["datetime"].dt.year <= 2019]
df["year"] = df["datetime"].dt.year
print(f"  Records after filtering (2000-2019): {len(df):,}")

def annual_maxima(df_in, dur_steps):
    rolled = df_in["depth_mm"].rolling(dur_steps).sum()
    df_in  = df_in.copy()
    df_in["rolled"] = rolled * (2 / dur_steps)
    return df_in.groupby("year")["rolled"].max().dropna()

annual_max = {}
for dur_min, label in zip(DURATIONS_MIN, DURATIONS_LABEL):
    annual_max[dur_min] = annual_maxima(df, dur_min // 30)

# ── Bootstrap ────────────────────────────────────────────────────────────────
print(f"Running {N_BOOT} bootstrap replicates ...")
rng = np.random.default_rng(SEED)
rows = []

for dur_min, label in zip(DURATIONS_MIN, DURATIONS_LABEL):
    am = annual_max[dur_min].values
    n  = len(am)
    boot_alpha = []

    for _ in range(N_BOOT):
        sample = rng.choice(am, size=n, replace=True)
        loc, scale = lmoments_gumbel(sample)
        ratios = []
        for T in RETURN_PERIODS:
            imerg_I = gumbel_quantile(T, loc, scale)
            if imerg_I > 0:
                ratios.append(OFFICIAL[dur_min][T] / imerg_I)
        if ratios:
            boot_alpha.append(np.mean(ratios))

    boot_alpha = np.array(boot_alpha)
    rows.append({
        "duration":      label,
        "alpha_mean":    round(float(np.mean(boot_alpha)), 2),
        "ci_lower_95":   round(float(np.percentile(boot_alpha, 2.5)), 2),
        "ci_upper_95":   round(float(np.percentile(boot_alpha, 97.5)), 2),
    })

out_df = pd.DataFrame(rows)
out_df.to_csv(OUT_PATH, index=False)
print(f"\nSaved -> {OUT_PATH}")
print(out_df.to_string(index=False))
