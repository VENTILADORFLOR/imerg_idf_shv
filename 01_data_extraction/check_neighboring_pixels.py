"""
check_neighboring_pixels.py

Check the IMERG pixel covering the study point and its 8 neighbouring pixels.
For each pixel, samples mean precipitation and flags whether the pixel
centroid is over water (using the JRC Global Surface Water layer).

Run this in your local environment (same setup as imerg_2000_2019.py).
Output: printed table to console
"""

import sys
from pathlib import Path
import ee

sys.path.append(str(Path(__file__).resolve().parents[1]))
from config.gee_config import GEE_PROJECT_ID, STUDY_POINT

# ── GEE initialisation ────────────────────────────────────────────────────────
ee.Initialize(project=GEE_PROJECT_ID)

lon, lat = STUDY_POINT
res = 0.1  # IMERG grid spacing in degrees

offsets = [
    (-1,  1), (0,  1), (1,  1),
    (-1,  0), (0,  0), (1,  0),
    (-1, -1), (0, -1), (1, -1),
]

imerg_mean = (
    ee.ImageCollection("NASA/GPM_L3/IMERG_V07")
    .filterDate("2000-06-01", "2001-01-01")   # one year only — faster
    .select("precipitation")
    .mean()
)

water = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence")
combined = imerg_mean.rename("precip").addBands(water)

# Build a FeatureCollection of all 9 points — one getInfo() call
features = []
for dlon, dlat in offsets:
    px_lon = lon + dlon * res
    px_lat = lat + dlat * res
    label  = f"({dlon:+d},{dlat:+d})" if not (dlon == 0 and dlat == 0) else "(0,0)STUDY"
    pt = ee.Geometry.Point([px_lon, px_lat])
    feat = ee.Feature(pt, {"label": label, "px_lon": px_lon, "px_lat": px_lat})
    features.append(feat)

fc = ee.FeatureCollection(features)

sampled = combined.sampleRegions(
    collection=fc,
    scale=5000,
    geometries=False,
)

print("Querying GEE (single request) ...")
result = sampled.getInfo()

print(f"\n{'Label':>12}  {'Lon':>10}  {'Lat':>9}  {'Mean precip':>13}  {'Water occ (%)':>15}")
print("-" * 68)

for feat in result["features"]:
    p = feat["properties"]
    label   = p.get("label", "?")
    px_lon  = p.get("px_lon", 0)
    px_lat  = p.get("px_lat", 0)
    precip  = p.get("precip")
    w_occ   = p.get("occurrence")
    p_str   = f"{precip:.4f}" if precip is not None else "N/A"
    w_str   = f"{w_occ:.1f}"  if w_occ   is not None else "N/A"
    print(f"{label:>12}  {px_lon:>10.4f}  {px_lat:>9.4f}  {p_str:>13}  {w_str:>15}")

print("\nWater occurrence > 50% = predominantly water pixel.")
