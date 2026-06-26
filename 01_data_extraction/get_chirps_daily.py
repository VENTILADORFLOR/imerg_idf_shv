"""
01_data_extraction/get_chirps_daily.py

Extract CHIRPS v3 daily precipitation for the Sihanoukville study point,
full available period 1981-01-01 to 2026-01-01, and export to Google Drive.

CHIRPS v3 (Climate Hazards Center InfraRed Precipitation with Station data):
    - Spatial resolution : 0.05 degrees (~5 km)
    - Temporal resolution: daily
    - Period             : 1981 to near-present
    - GEE collection     : UCSB-CHC/CHIRPS/V3/DAILY_RNL
    - Reference          : Funk et al. (2015), Scientific Data 2:150066

Note on collection choice:
    Two daily CHIRPS v3 products are available on GEE:
        DAILY_SAT  uses NASA IMERG Late V07 to partition pentadal totals
                   into daily amounts. IMERG starts in 2000, so DAILY_SAT
                   is unavailable before ~1998 despite CHIRPS itself
                   starting in 1981.
        DAILY_RNL  uses ECMWF ERA5 reanalysis to partition pentadal totals
                   into daily amounts. ERA5 covers 1940 to present, so
                   DAILY_RNL provides the full 1981-to-present record.
    DAILY_RNL is used here to obtain the complete 1981-2025 time series
    required for the 45-year Gumbel comparison in the manuscript.

Output (Google Drive -> data/raw/):
    rainfall_point_chirps_1981_2025.csv

Columns:
    date                : YYYY-MM-DD
    chirps_precipitation: daily rainfall (mm/day)

Note:
    CHIRPS is daily only — no sub-daily resolution.
    This extract is used exclusively for 24-hour Gumbel comparison
    against IMERG and the official IDF benchmark.
"""

import sys
from pathlib import Path
import ee
import time

sys.path.append(str(Path(__file__).resolve().parents[1]))
from config.gee_config import GEE_PROJECT_ID, STUDY_POINT

# ── GEE initialisation ────────────────────────────────────────────────────────
ee.Initialize(project=GEE_PROJECT_ID)

# ── Study point ───────────────────────────────────────────────────────────────
point = ee.Geometry.Point(STUDY_POINT)

# ── CHIRPS v3 DAILY_RNL collection ───────────────────────────────────────────
# ERA5-based daily partitioning; covers full 1981-to-present record.
chirps = (
    ee.ImageCollection("UCSB-CHC/CHIRPS/V3/DAILY_RNL")
    .filterDate("1981-01-01", "2026-01-01")
    .select("precipitation")
)

# ── Feature mapping ───────────────────────────────────────────────────────────
def image_to_feature(image):
    value = image.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=point,
        scale=5000,
        maxPixels=1e9,
    ).get("precipitation")
    return ee.Feature(None, {
        "date":                 image.date().format("YYYY-MM-dd"),
        "chirps_precipitation": value,
    })

features = ee.FeatureCollection(chirps.map(image_to_feature))

# ── Export to Google Drive ────────────────────────────────────────────────────
task = ee.batch.Export.table.toDrive(
    collection=features,
    description="Rainfall_Point_CHIRPS_1981_2025",
    folder="GEE_Exports",
    fileNamePrefix="rainfall_point_chirps_1981_2025",
    fileFormat="CSV",
)

task.start()
print("CHIRPS v3 (DAILY_RNL) export task submitted.")

while True:
    status = task.status()
    state = status["state"]
    if state == "COMPLETED":
        print("Export completed. Download from Google Drive -> data/raw/")
        break
    elif state == "FAILED":
        print(f"Export failed: {status['error_message']}")
        break
    else:
        print(f"Task state: {state} ...")
        time.sleep(30)
