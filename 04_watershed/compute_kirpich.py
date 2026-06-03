"""
04_watershed/compute_kirpich.py

Compute the main flow path length and average slope within the delineated
catchment, then estimate the time of concentration using the Kirpich (1940)
formula. Results are printed to console and saved as CSV.

Inputs:
    data/watershed/catchment.shp
    data/watershed/flow_paths.shp
    data/watershed/dem_utm48n.tif   (fallback: data/raw/sihanoukville_alos30m_dem.tif)

Output:
    data/processed/kirpich_tc.csv
"""

from pathlib import Path
import numpy as np
import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.mask import mask as rio_mask

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parents[1]
CATCHMENT  = ROOT / "data" / "watershed" / "catchment.shp"
FLOW_PATHS = ROOT / "data" / "watershed" / "flow_paths.shp"
DEM_UTM    = ROOT / "data" / "watershed" / "dem_utm48n.tif"
DEM_RAW    = ROOT / "data" / "raw"       / "sihanoukville_alos30m_dem.tif"
OUT_PATH   = ROOT / "data" / "processed" / "kirpich_tc.csv"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

dem_path = DEM_UTM if DEM_UTM.exists() else DEM_RAW
print(f"Using DEM: {dem_path}")

# ── Flow path length ──────────────────────────────────────────────────────────
catch_gdf    = gpd.read_file(CATCHMENT)
flow_gdf     = gpd.read_file(FLOW_PATHS)
flow_clipped = gpd.clip(flow_gdf, catch_gdf)
flow_clipped["length_m"] = flow_clipped.geometry.length.round(1)
L_m = float(flow_clipped["length_m"].sum())
print(f"Total flow path length: {L_m:.1f} m")

# ── Elevation difference ──────────────────────────────────────────────────────
catchment_geom = [catch_gdf.geometry.iloc[0].__geo_interface__]
with rasterio.open(dem_path) as src:
    src_nodata = src.nodata if src.nodata is not None else -9999
    dem_data, _ = rio_mask(src, catchment_geom, crop=True, nodata=src_nodata)
    dem_arr = dem_data[0].astype(float)
    dem_arr[dem_arr == src_nodata] = np.nan

elev_max = float(np.nanmax(dem_arr))
elev_min = float(np.nanmin(dem_arr))
delta_H  = elev_max - elev_min
S        = delta_H / L_m
print(f"Elevation: max={elev_max:.1f} m, min={elev_min:.1f} m, ΔH={delta_H:.1f} m")
print(f"Average slope S = {S:.4f} m/m ({S*100:.2f}%)")

# ── Kirpich (1940) ────────────────────────────────────────────────────────────
Tc_min = 0.0195 * (L_m ** 0.77) * (S ** -0.385)
print(f"Kirpich Tc = {Tc_min:.1f} minutes")

# ── Save CSV ──────────────────────────────────────────────────────────────────
out_df = pd.DataFrame([{
    "flow_path_length_m":   round(L_m, 1),
    "elev_max_m":           round(elev_max, 1),
    "elev_min_m":           round(elev_min, 1),
    "delta_H_m":            round(delta_H, 1),
    "slope_m_per_m":        round(S, 4),
    "slope_pct":            round(S * 100, 2),
    "kirpich_tc_min":       round(Tc_min, 1),
}])
out_df.to_csv(OUT_PATH, index=False)
print(f"Saved -> {OUT_PATH}")
print(out_df.to_string(index=False))
