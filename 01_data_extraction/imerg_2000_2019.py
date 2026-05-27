"""
01_data_extraction/imerg_2000_2019.py

Extract GPM IMERG V07 half-hourly precipitation for the Sihanoukville study point,
period 2000-06-01 to 2020-01-01, and export to Google Drive as CSV.

Output (Google Drive → data/raw/):
    rainfall_point_imerg_2000_2019.csv

Columns:
    datetime            : YYYY-MM-DD HH:mm (UTC)
    gpm_precipitation   : precipitation rate (mm/hr)
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

# ── IMERG collection ──────────────────────────────────────────────────────────
imerg = (
    ee.ImageCollection("NASA/GPM_L3/IMERG_V07")
    .filterDate("2000-06-01", "2020-01-01")
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
        "datetime": image.date().format("YYYY-MM-dd HH:mm"),
        "gpm_precipitation": value,
    })

features = ee.FeatureCollection(imerg.map(image_to_feature))

# ── Export to Google Drive ────────────────────────────────────────────────────
task = ee.batch.Export.table.toDrive(
    collection=features,
    description="Rainfall_Point_IMERG_2000_2019",
    folder="GEE_Exports",
    fileNamePrefix="rainfall_point_imerg_2000_2019",
    fileFormat="CSV",
)

task.start()
print("Export task submitted.")

while True:
    status = task.status()
    state = status["state"]
    if state == "COMPLETED":
        print("Export completed. Download from Google Drive → data/raw/")
        break
    elif state == "FAILED":
        print(f"Export failed: {status['error_message']}")
        break
    else:
        print(f"Task state: {state} ...")
        time.sleep(30)
