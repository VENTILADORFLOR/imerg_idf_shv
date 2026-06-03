"""
04_watershed/delineate_watershed.py

Automatically delineate catchment boundary and flow path network
from ALOS 30m DEM, using the drainage node as the pour point.

CRS strategy:
    ALL inputs reprojected to EPSG:32648 (WGS 1984 UTM Zone 48N) at load time.
    ALL outputs saved in EPSG:32648.

Input:
    data/raw/sihanoukville_alos30m_dem.tif
    data/drainage_node/drainage_node.shp

Output (data/watershed/):
    catchment.shp     [EPSG:32648]
    flow_paths.shp    [EPSG:32648]

Requirements:
    pip install pysheds geopandas fiona shapely rasterio
"""

from pathlib import Path
import numpy as np
import geopandas as gpd
from shapely.geometry import shape, LineString
from shapely.ops import unary_union
from pysheds.grid import Grid
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import warnings
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT   = Path(__file__).resolve().parents[1]
DEM    = ROOT / "data" / "raw"           / "sihanoukville_alos30m_dem.tif"
NODE   = ROOT / "data" / "drainage_node" / "drainage_node.shp"
OUTPUT = ROOT / "data" / "watershed"
OUTPUT.mkdir(parents=True, exist_ok=True)

TARGET_CRS       = "EPSG:32648"   # WGS 1984 UTM Zone 48N — metric, covers Cambodia
DEM_UTM          = OUTPUT / "dem_utm48n.tif"
STREAM_THRESHOLD = 500            # cells; lower = more detail, higher = main channels only

# ── Helper: reproject any GeoDataFrame ───────────────────────────────────────
def load_shp(path, label=""):
    gdf = gpd.read_file(path)
    src_crs = gdf.crs
    if str(src_crs) != TARGET_CRS:
        gdf = gdf.to_crs(TARGET_CRS)
        print(f"  {label}: {src_crs} → {TARGET_CRS}")
    else:
        print(f"  {label}: already {TARGET_CRS}")
    return gdf

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Reproject DEM to UTM 48N
# ══════════════════════════════════════════════════════════════════════════════
print("Step 1: Reprojecting DEM to EPSG:32648 ...")
with rasterio.open(DEM) as src:
    if src.crs.to_epsg() == 32648:
        print("  DEM already in EPSG:32648, skipping.")
        DEM_UTM = DEM
    else:
        transform, width, height = calculate_default_transform(
            src.crs, TARGET_CRS, src.width, src.height, *src.bounds)
        meta = src.meta.copy()
        meta.update({"crs": TARGET_CRS, "transform": transform,
                     "width": width, "height": height})
        with rasterio.open(DEM_UTM, "w", **meta) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=TARGET_CRS,
                    resampling=Resampling.bilinear,
                )
        print(f"  Saved reprojected DEM → {DEM_UTM}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Load and reproject drainage node
# ══════════════════════════════════════════════════════════════════════════════
print("Step 2: Loading inputs (auto-reproject to EPSG:32648) ...")
node_gdf = load_shp(NODE, "drainage_node")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Condition DEM
# ══════════════════════════════════════════════════════════════════════════════
print("Step 3: Conditioning DEM ...")
grid = Grid.from_raster(str(DEM_UTM))
dem  = grid.read_raster(str(DEM_UTM))

pit_filled = grid.fill_pits(dem)
flooded    = grid.fill_depressions(pit_filled)
inflated   = grid.resolve_flats(flooded)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Flow direction and accumulation
# ══════════════════════════════════════════════════════════════════════════════
print("Step 4: Computing flow direction and accumulation ...")
fdir = grid.flowdir(inflated)
acc  = grid.accumulation(fdir)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Snap pour point to stream
# ══════════════════════════════════════════════════════════════════════════════
print("Step 5: Snapping pour point ...")
pt = node_gdf.geometry.iloc[0]
x_snap, y_snap = grid.snap_to_mask(acc > STREAM_THRESHOLD, (pt.x, pt.y))
offset = ((x_snap - pt.x)**2 + (y_snap - pt.y)**2) ** 0.5
print(f"  Original : ({pt.x:.1f}, {pt.y:.1f}) m")
print(f"  Snapped  : ({x_snap:.1f}, {y_snap:.1f}) m  (offset {offset:.1f} m)")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Delineate catchment
# ══════════════════════════════════════════════════════════════════════════════
print("Step 6: Delineating catchment ...")
catch = grid.catchment(x=x_snap, y=y_snap, fdir=fdir, xytype='coordinate')
grid.clip_to(catch)
catch_view = grid.view(catch).astype(np.uint8)   # rasterio polygonize requires explicit dtype
polys = [shape(s) for s, v in grid.polygonize(catch_view) if v == 1]

if not polys:
    raise RuntimeError(
        "No catchment polygon generated.\n"
        "Check that drainage_node.shp falls within the DEM extent.")

catch_gdf = gpd.GeoDataFrame(geometry=[unary_union(polys)], crs=TARGET_CRS)
area_km2  = catch_gdf.geometry.area.iloc[0] / 1e6
catch_gdf["area_km2"] = round(area_km2, 3)
catch_gdf.to_file(OUTPUT / "catchment.shp")
print(f"  Catchment area : {area_km2:.3f} km²  (paper reports 3.43 km²)")
print(f"  Saved → {OUTPUT / 'catchment.shp'}  [EPSG:32648]")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Extract stream network
# ══════════════════════════════════════════════════════════════════════════════
print("Step 7: Extracting stream network ...")
streams = grid.extract_river_network(fdir=fdir, mask=acc > STREAM_THRESHOLD)
lines   = [LineString(b['geometry']['coordinates'])
           for b in streams['features']
           if len(b['geometry']['coordinates']) >= 2]

if not lines:
    raise RuntimeError("No streams extracted. Try lowering STREAM_THRESHOLD.")

flow_gdf = gpd.GeoDataFrame(geometry=lines, crs=TARGET_CRS)
flow_gdf["length_m"] = flow_gdf.geometry.length.round(1)
flow_gdf.to_file(OUTPUT / "flow_paths.shp")
print(f"  Stream segments : {len(flow_gdf)}")
print(f"  Saved → {OUTPUT / 'flow_paths.shp'}  [EPSG:32648]")

print("\n✓ Done. All outputs in EPSG:32648:")
print(f"  {OUTPUT / 'catchment.shp'}")
print(f"  {OUTPUT / 'flow_paths.shp'}")
