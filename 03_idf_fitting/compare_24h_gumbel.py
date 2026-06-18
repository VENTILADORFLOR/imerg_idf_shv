"""
03_idf_fitting/compare_24h_gumbel.py

Compare 24-hour Gumbel frequency curves from five sources:

    1. Official IDF benchmark (Department of Hydrology, Cambodia, 1999-2019, n=20)
    2. CHIRPS v3 (1999-2019, n=21) — same window as official, isolates 1999 outlier effect
    3. IMERG V07 (2000-2019, n=20) — primary analysis period, IMERG starts 2000-06
    4. CHIRPS v3 (1981-2025, n=44) — long record, tests whether official 24h depths are elevated
    5. IMERG V07 (2000-2024, n=25) — extended record for temporal stability check

Two analytical purposes:
    Purpose 1 (same-window comparison, lines 1-3): assess satellite vs official bias
              over an equivalent period; CHIRPS 1999-2019 bridges the one-year gap
              between IMERG (starts 2000) and the official record (starts 1999).
    Purpose 2 (long-record comparison, lines 4-5): assess whether the official 20-year
              record — which includes the 1999 extreme event — systematically overestimates
              24-hour design depths relative to longer satellite records.

Inputs:
    data/raw/rainfall_point_imerg_2000_2019.csv   (existing)
    data/raw/rainfall_point_imerg_2020_2024.csv   (existing)
    data/raw/rainfall_point_chirps_1981_2025.csv  (from get_chirps_daily.py)

Outputs:
    data/figures/gumbel_24h_comparison.png
    data/processed/gumbel_24h_comparison.csv

References:
    Funk et al. (2015). The climate hazards infrared precipitation with
    stations. Scientific Data 2:150066.
    CHIRPS v3: UCSB-CHC/CHIRPS/V3/DAILY_RNL (GEE collection).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import gumbel_r
from scipy.optimize import curve_fit

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT    = Path(__file__).resolve().parents[1]
RAW     = ROOT / "data" / "raw"
OUT_CSV = ROOT / "data" / "processed" / "gumbel_24h_comparison.csv"
OUT_FIG = ROOT / "data" / "figures"   / "gumbel_24h_comparison.png"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
OUT_FIG.parent.mkdir(parents=True, exist_ok=True)

# ── Official 24h design depths (mm) ──────────────────────────────────────────
# Tabulated values from official IDF benchmark; confirmed to lie perfectly on
# a single Gumbel curve (reverse-engineered from Formula A).
OFFICIAL_P24H = {
    2.33: 129.9,
    5:    170.2,
    10:   203.0,
    25:   244.5,
    50:   275.2,
    100:  305.8,
}

# ── Return periods ────────────────────────────────────────────────────────────
T_PLOT    = np.logspace(np.log10(1.5), np.log10(200), 200)
T_MARKERS = [2.33, 5, 10, 25, 50, 100]

# ── Helpers ───────────────────────────────────────────────────────────────────
def fit_gumbel(annual_max):
    """Fit Gumbel by MLE; return (quantile_fn, loc, scale)."""
    loc, scale = gumbel_r.fit(annual_max)
    def quantile(T):
        return gumbel_r.ppf(1 - 1 / np.asarray(T), loc=loc, scale=scale)
    return quantile, loc, scale

def gumbel_quantile_model(T, loc, scale):
    return loc - scale * np.log(-np.log(1 - 1 / T))

def imerg_rolling_24h(df_raw, year_start, year_end):
    """
    From merged IMERG half-hourly DataFrame, compute 24-hour rolling annual
    maxima for a specified year range [year_start, year_end] inclusive.
    """
    df = df_raw.copy()
    df = df.sort_values("datetime").set_index("datetime")
    df = df.resample("30min").mean(numeric_only=True)
    df["depth_mm_30min"] = df["gpm_precipitation"] * 0.5
    df["rolling_24h"] = df["depth_mm_30min"].rolling(48).sum()
    df["year"] = df.index.year
    df = df[(df["year"] >= year_start) & (df["year"] <= year_end)]
    am = df.groupby("year")["rolling_24h"].max().dropna()
    return am

def chirps_annual_max(df_raw, year_start, year_end):
    """Extract CHIRPS annual maxima for a specified year range."""
    df = df_raw.copy()
    df = df[(df["year"] >= year_start) & (df["year"] <= year_end)]
    am = df.groupby("year")["chirps_precipitation"].max().dropna()
    return am

# ══════════════════════════════════════════════════════════════════════════════
# Load raw data
# ══════════════════════════════════════════════════════════════════════════════
print("Loading IMERG ...")
df_a = pd.read_csv(RAW / "rainfall_point_imerg_2000_2019.csv")
df_b = pd.read_csv(RAW / "rainfall_point_imerg_2020_2024.csv")
df_imerg_raw = pd.concat([df_a, df_b], ignore_index=True)
df_imerg_raw["datetime"] = pd.to_datetime(df_imerg_raw["datetime"])
df_imerg_raw = df_imerg_raw.dropna(subset=["gpm_precipitation"])
df_imerg_raw = df_imerg_raw.drop_duplicates(subset=["datetime"])

print("Loading CHIRPS ...")
df_chirps_raw = pd.read_csv(RAW / "rainfall_point_chirps_1981_2025.csv")
df_chirps_raw["date"] = pd.to_datetime(df_chirps_raw["date"])
df_chirps_raw = df_chirps_raw.dropna(subset=["chirps_precipitation"])
df_chirps_raw = df_chirps_raw[df_chirps_raw["chirps_precipitation"] >= 0]
df_chirps_raw["year"] = df_chirps_raw["date"].dt.year

# ══════════════════════════════════════════════════════════════════════════════
# 1. Official — reverse-fit Gumbel from tabulated design depths
# ══════════════════════════════════════════════════════════════════════════════
print("\nFitting official Gumbel ...")
T_known = np.array(list(OFFICIAL_P24H.keys()))
P_known = np.array(list(OFFICIAL_P24H.values()))
popt, _ = curve_fit(gumbel_quantile_model, T_known, P_known, p0=[100, 50])
off_loc, off_scale = popt
def official_q(T):
    return gumbel_quantile_model(np.asarray(T), off_loc, off_scale)
print(f"  Official   Gumbel: loc={off_loc:.2f}, scale={off_scale:.2f}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. CHIRPS 1999-2019 (same window as official; includes 1999 extreme event)
# ══════════════════════════════════════════════════════════════════════════════
chirps_am_1999_2019 = chirps_annual_max(df_chirps_raw, 1999, 2019)
n_c19 = len(chirps_am_1999_2019)
chirps_q_1999_2019, c19_loc, c19_scale = fit_gumbel(chirps_am_1999_2019.values)
print(f"  CHIRPS 1999-2019   Gumbel: loc={c19_loc:.2f}, scale={c19_scale:.2f}  (n={n_c19})")

# ══════════════════════════════════════════════════════════════════════════════
# 3. IMERG 2000-2019 (primary analysis period)
# ══════════════════════════════════════════════════════════════════════════════
imerg_am_2000_2019 = imerg_rolling_24h(df_imerg_raw, 2000, 2019)
n_i19 = len(imerg_am_2000_2019)
imerg_q_2000_2019, i19_loc, i19_scale = fit_gumbel(imerg_am_2000_2019.values)
print(f"  IMERG 2000-2019    Gumbel: loc={i19_loc:.2f}, scale={i19_scale:.2f}  (n={n_i19})")

# ══════════════════════════════════════════════════════════════════════════════
# 4. CHIRPS 1981-2025 (long record)
# ══════════════════════════════════════════════════════════════════════════════
chirps_am_full = chirps_annual_max(df_chirps_raw, 1981, 2025)
n_cf = len(chirps_am_full)
chirps_q_full, cf_loc, cf_scale = fit_gumbel(chirps_am_full.values)
print(f"  CHIRPS 1981-2025   Gumbel: loc={cf_loc:.2f}, scale={cf_scale:.2f}  (n={n_cf})")

# ══════════════════════════════════════════════════════════════════════════════
# 5. IMERG 2000-2024 (extended record)
# ══════════════════════════════════════════════════════════════════════════════
imerg_am_full = imerg_rolling_24h(df_imerg_raw, 2000, 2024)
n_if = len(imerg_am_full)
imerg_q_full, if_loc, if_scale = fit_gumbel(imerg_am_full.values)
print(f"  IMERG 2000-2024    Gumbel: loc={if_loc:.2f}, scale={if_scale:.2f}  (n={n_if})")

# ══════════════════════════════════════════════════════════════════════════════
# Export comparison table
# ══════════════════════════════════════════════════════════════════════════════
rows = []
for T in T_MARKERS:
    rows.append({
        "return_period_year":          T,
        "official_P24h_mm":            round(official_q(T),                   1),
        "chirps_1999_2019_P24h_mm":    round(float(chirps_q_1999_2019(T)),    1),
        "imerg_2000_2019_P24h_mm":     round(float(imerg_q_2000_2019(T)),     1),
        "chirps_1981_2025_P24h_mm":    round(float(chirps_q_full(T)),         1),
        "imerg_2000_2024_P24h_mm":     round(float(imerg_q_full(T)),          1),
    })

out_df = pd.DataFrame(rows)
out_df.to_csv(OUT_CSV, index=False)
print(f"\nSaved -> {OUT_CSV}")
print(out_df.to_string(index=False))

# ══════════════════════════════════════════════════════════════════════════════
# Plot
# ══════════════════════════════════════════════════════════════════════════════
print("\nPlotting ...")
fig, ax = plt.subplots(figsize=(9, 5.5))

# ── Five Gumbel curves ────────────────────────────────────────────────────────
# Official
ax.plot(T_PLOT, official_q(T_PLOT),
        color="#D62728", lw=2.2, ls="-",
        label=f"Official benchmark (1999–2019, n=21)")

# CHIRPS 1999-2019 (same window, dashed green)
ax.plot(T_PLOT, chirps_q_1999_2019(T_PLOT),
        color="#2CA02C", lw=1.6, ls="--",
        label=f"CHIRPS v3 (1999–2019, n={n_c19})")

# IMERG 2000-2019 (primary period, dashed blue)
ax.plot(T_PLOT, imerg_q_2000_2019(T_PLOT),
        color="#1F77B4", lw=1.6, ls="--",
        label=f"GPM IMERG V07 (2000–2019, n={n_i19})")

# CHIRPS 1981-2025 (long record, solid green)
ax.plot(T_PLOT, chirps_q_full(T_PLOT),
        color="#2CA02C", lw=2.2, ls="-",
        label=f"CHIRPS v3 (1981–2025, n={n_cf})")

# IMERG 2000-2024 (extended, solid blue)
ax.plot(T_PLOT, imerg_q_full(T_PLOT),
        color="#1F77B4", lw=2.2, ls="-",
        label=f"GPM IMERG V07 (2000–2024, n={n_if})")

# ── Scatter: official tabulated points ───────────────────────────────────────
ax.scatter(list(OFFICIAL_P24H.keys()), list(OFFICIAL_P24H.values()),
           color="#D62728", s=45, zorder=6)

# ── Scatter: IMERG 2000-2019 annual maxima (Weibull positions) ────────────────
T_w = (n_i19 + 1) / np.arange(1, n_i19 + 1)
ax.scatter(sorted(T_w), sorted(imerg_am_2000_2019.values),
           color="#1F77B4", s=18, alpha=0.45, zorder=4, marker="^")

# ── Scatter: CHIRPS 1999-2019 annual maxima ───────────────────────────────────
T_wc = (n_c19 + 1) / np.arange(1, n_c19 + 1)
ax.scatter(sorted(T_wc), sorted(chirps_am_1999_2019.values),
           color="#2CA02C", s=18, alpha=0.45, zorder=4, marker="s")

# ── Reference lines ───────────────────────────────────────────────────────────
for T_ref in [5, 10, 25]:
    ax.axvline(T_ref, color="grey", lw=0.6, ls="--", alpha=0.5)

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xscale("log")
ax.set_xlabel("Return period (years)", fontsize=11)
ax.set_ylabel("24-hour rainfall depth (mm)", fontsize=11)
ax.set_title(
    "24-hour Gumbel frequency curves: official benchmark, IMERG, and CHIRPS\n"
    "Solid lines = full record; dashed lines = 1999/2000–2019 window",
    fontsize=10,
)
ax.legend(fontsize=8.5, framealpha=0.9, loc="upper left")
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.set_xticks([2, 5, 10, 25, 50, 100])
ax.grid(True, which="both", alpha=0.3)
ax.set_xlim(1.5, 200)

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved -> {OUT_FIG}")
