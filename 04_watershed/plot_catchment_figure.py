"""
04_watershed/plot_catchment_figure.py

Reproduce the catchment analysis figure from the paper.
ALL input layers are reprojected to EPSG:32648 at load time,
regardless of their original CRS (WGS84 or UTM 48N).

Input:
    data/raw/sihanoukville_alos30m_dem.tif      ← hillshade base
    data/study_area/study_area.shp              ← red site boundary
    data/depression/depression.shp              ← blue depression areas
    data/drainage_node/drainage_node.shp        ← orange pour point
    data/main_road/main_road.shp                ← cyan main road
    data/roads/gis_osm_roads_free_1.shp         ← light blue road grid
    data/watershed/catchment.shp                ← green catchment area
    data/watershed/flow_paths.shp               ← magenta flow paths

Output:
    data/figures/catchment_figure.png

Requirements:
    pip install matplotlib geopandas rasterio numpy
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import geopandas as gpd
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.enums import Resampling as RResampling
import warnings
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parents[1]
DATA        = ROOT / "data"
FIGURES     = DATA / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

TARGET_CRS = "EPSG:32648"   # WGS 1984 UTM Zone 48N

# ── Helper: load and reproject any shapefile ──────────────────────────────────
def load_shp(path, label=""):
    gdf = gpd.read_file(path)
    if str(gdf.crs) != TARGET_CRS:
        print(f"  Reprojecting {label}: {gdf.crs} → {TARGET_CRS}")
        gdf = gdf.to_crs(TARGET_CRS)
    else:
        print(f"  {label}: already {TARGET_CRS}")
    return gdf

# ── Helper: reproject raster and return hillshade array ──────────────────────
def make_hillshade(dem_path):
    tmp = FIGURES / "dem_utm_tmp.tif"
    with rasterio.open(dem_path) as src:
        if src.crs.to_epsg() == 32648:
            with rasterio.open(dem_path) as f:
                elev = f.read(1).astype(float)
                transform = f.transform
                bounds = f.bounds
        else:
            t, w, h = calculate_default_transform(
                src.crs, TARGET_CRS, src.width, src.height, *src.bounds)
            meta = src.meta.copy()
            meta.update({"crs": TARGET_CRS, "transform": t, "width": w, "height": h})
            with rasterio.open(tmp, "w", **meta) as dst:
                reproject(source=rasterio.band(src, 1),
                          destination=rasterio.band(dst, 1),
                          src_transform=src.transform, src_crs=src.crs,
                          dst_transform=t, dst_crs=TARGET_CRS,
                          resampling=Resampling.bilinear)
            with rasterio.open(tmp) as f:
                elev = f.read(1).astype(float)
                transform = f.transform
                bounds = f.bounds

    # ── Classic single-source hillshade (315° NW, altitude 45°) ─────────────
    # Multi-source was causing brightness overload on flat areas.
    # Single NW source gives the natural light/shadow gradient humans expect.
    res      = transform.a          # pixel size in metres
    altitude = 45 * np.pi / 180
    azimuth  = 315 * np.pi / 180

    # Replace nodata / negative fill values with local median before gradients
    elev = np.where(elev < -9000, np.nan, elev)
    elev_filled = np.where(np.isnan(elev), np.nanmedian(elev), elev)

    dzdx = np.gradient(elev_filled, res, axis=1)
    dzdy = np.gradient(elev_filled, res, axis=0)
    slope  = np.arctan(np.sqrt(dzdx**2 + dzdy**2))
    aspect = np.arctan2(-dzdy, dzdx)

    hs = (np.sin(altitude) * np.cos(slope)
          + np.cos(altitude) * np.sin(slope) * np.cos(azimuth - aspect))
    hs = np.clip(hs, 0, 1)

    # ── Linear mid-range mapping — locks contrast to creamy grey zone ────────
    # Deepest shadow → 0.35 (light grey), brightest plain → 0.95 (off-white).
    # Eliminates dead-black patches without any percentile clipping.
    hs_min = hs.min()
    hs_max = hs.max()
    hs = 0.35 + (hs - hs_min) * (0.95 - 0.35) / (hs_max - hs_min + 1e-9)

    return hs, bounds

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load all layers (auto-reproject to UTM 48N)
# ══════════════════════════════════════════════════════════════════════════════
print("Loading layers ...")
study_area   = load_shp(DATA / "study_area"    / "study_area.shp",              "study_area")
depression   = load_shp(DATA / "depression"    / "depression.shp",              "depression")
node         = load_shp(DATA / "drainage_node" / "drainage_node.shp",           "drainage_node")
main_road    = load_shp(DATA / "main_road"     / "main_road.shp",               "main_road")
roads        = load_shp(DATA / "roads"         / "gis_osm_roads_free_1.shp",    "roads")
catchment      = load_shp(DATA / "watershed"      / "catchment.shp",              "catchment")
flow_paths     = load_shp(DATA / "watershed"      / "flow_paths.shp",             "flow_paths")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Hillshade from DEM
# ══════════════════════════════════════════════════════════════════════════════
print("Generating hillshade ...")
hillshade, dem_bounds = make_hillshade(DATA / "raw" / "sihanoukville_alos30m_dem.tif")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Set map extent — union of study_area + catchment, A4 landscape aspect ratio
# ══════════════════════════════════════════════════════════════════════════════
from shapely.geometry import box  # noqa: already imported above via shapely

# Combine study_area and catchment bounds to ensure both are fully visible
combined_bounds = np.array([
    min(study_area.total_bounds[0], catchment.total_bounds[0]),
    min(study_area.total_bounds[1], catchment.total_bounds[1]),
    max(study_area.total_bounds[2], catchment.total_bounds[2]),
    max(study_area.total_bounds[3], catchment.total_bounds[3]),
])
buf  = 600                             # metres padding around combined extent
cx0  = combined_bounds[0] - buf
cy0  = combined_bounds[1] - buf
cx1  = combined_bounds[2] + buf
cy1  = combined_bounds[3] + buf

# Force A4 landscape aspect ratio (297 / 210 ≈ 1.4143)
A4_RATIO = 297 / 210
data_w   = cx1 - cx0
data_h   = cy1 - cy0
if data_w / data_h < A4_RATIO:        # too tall → widen
    extra = data_h * A4_RATIO - data_w
    cx0  -= extra / 2
    cx1  += extra / 2
else:                                  # too wide → heighten
    extra = data_w / A4_RATIO - data_h
    cy0  -= extra / 2
    cy1  += extra / 2

xmin, ymin, xmax, ymax = cx0, cy0, cx1, cy1

# ══════════════════════════════════════════════════════════════════════════════
# 4. Plot
# ══════════════════════════════════════════════════════════════════════════════
print("Plotting ...")
fig, ax = plt.subplots(figsize=(16.54, 11.69), dpi=300)  # A4 landscape @ 300 dpi (print quality)

# Hillshade extent in UTM metres
hs_extent = [dem_bounds.left, dem_bounds.right, dem_bounds.bottom, dem_bounds.top]
ax.imshow(hillshade, extent=hs_extent, cmap="gray", vmin=0, vmax=1, alpha=0.9,
          origin="upper", zorder=0, interpolation="bilinear")

# Catchment area (green semi-transparent)
catchment.plot(ax=ax, facecolor="#90EE90", edgecolor="#228B22",
               alpha=0.45, linewidth=1.2, zorder=1)

# Roads — gis_osm_roads_free_1.shp (Geofabrik / OpenStreetMap) (light cyan road grid)
roads.plot(ax=ax, color="#00CED1", linewidth=0.8, alpha=0.8, zorder=2)

# Flow paths (magenta)
flow_paths.plot(ax=ax, color="#FF00FF", linewidth=0.8, zorder=3)

# Main road (cyan, bold, on top of roads)
main_road.plot(ax=ax, color="#00CED1", linewidth=2.0, zorder=4)

# Depression (blue fill)
depression.plot(ax=ax, facecolor="#6495ED", edgecolor="#4169E1",
                alpha=0.7, linewidth=0.8, zorder=5)

# Study area boundary (red bold)
study_area.plot(ax=ax, facecolor="none", edgecolor="red",
                linewidth=2.5, zorder=6)

# Drainage node (orange circle)
node.plot(ax=ax, color="#FFA500", markersize=80, zorder=7)

# ── Labels ────────────────────────────────────────────────────────────────────
# Catchment area annotation
catch_centroid = catchment.geometry.centroid.iloc[0]
area_km2 = catchment["area_km2"].iloc[0] if "area_km2" in catchment.columns else 0.0
ax.text(catch_centroid.x - 400, 
        catch_centroid.y + 470, 
        f"External catchment = {area_km2:.2f} km²",
        ha="center", va="center", fontsize=9.5, color="#1a1a1a",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.5, ec="none"))

# Critical drainage node label
node_pt = node.geometry.iloc[0]
ax.annotate("Critical drainage node",
            xy=(node_pt.x, node_pt.y),
            xytext=(node_pt.x + 120, node_pt.y + 80),
            fontsize=9, color="black",
            arrowprops=dict(arrowstyle="-", color="black", lw=0.8))

# ── Map extent ────────────────────────────────────────────────────────────────
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect("equal")
ax.axis("off")

# ── North arrow ───────────────────────────────────────────────────────────────
arrow_x = xmax - (xmax - xmin) * 0.05
arrow_y = ymax - (ymax - ymin) * 0.08
ax.annotate("", xy=(arrow_x, arrow_y),
            xytext=(arrow_x, arrow_y - (ymax - ymin) * 0.05),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5))
ax.text(arrow_x, arrow_y + (ymax - ymin) * 0.01, "N",
        ha="center", va="bottom", fontsize=11, fontweight="bold")

# ── Scale bar (300 m) ─────────────────────────────────────────────────────────
sb_len  = 300          # metres
sb_x0   = xmax - (xmax - xmin) * 0.22
sb_y0   = ymin + (ymax - ymin) * 0.04
ax.plot([sb_x0, sb_x0 + sb_len], [sb_y0, sb_y0], "k-", lw=2)
ax.plot([sb_x0, sb_x0],          [sb_y0 - 20, sb_y0 + 20], "k-", lw=2)
ax.plot([sb_x0 + sb_len, sb_x0 + sb_len], [sb_y0 - 20, sb_y0 + 20], "k-", lw=2)
ax.text(sb_x0 + sb_len / 2, sb_y0 - 60, "300 m",
        ha="center", va="top", fontsize=8)
ax.text(sb_x0, sb_y0 - 60, "0", ha="center", va="top", fontsize=8)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_elements = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#FFA500",
           markeredgewidth=0, markersize=8, label="Drainage node"),
    Line2D([0], [0], color="#00CED1", lw=2,   label="Main road"),
    Line2D([0], [0], color="#00CED1", lw=0.8, alpha=0.8, label="Roads"),
    Line2D([0], [0], color="#FF00FF", lw=1,   label="Flow paths"),
    Line2D([0], [0], color="red",     lw=2.5, label="Site boundary"),
    mpatches.Patch(fc="#6495ED", ec="#4169E1", alpha=0.7, label="Depression"),
    mpatches.Patch(fc="#90EE90", ec="#228B22", alpha=0.6, label="Catchment area"),
]
ax.legend(handles=legend_elements, loc="lower left",
          fontsize=8, framealpha=0.85, frameon=True)

# ── Save ──────────────────────────────────────────────────────────────────────
out_path = FIGURES / "catchment_figure.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"\n✓ Saved → {out_path}")
