"""
01_data_extraction/get_alos_dem.py

Export ALOS AW3D30 V4.1 DSM for the Sihanoukville study region
from Google Earth Engine to Google Drive as GeoTIFF.

Output (Google Drive → data/raw/):
    sihanoukville_alos30m_dem.tif

Notes:
    - AW3D30 V4.1 (April 2024): corrects anomalies found in V3.1/V3.2 for
      19,051 tiles globally; band name remains DSM, native resolution 30 m.
    - scale=30 matches the native resolution; using 12.5 m would trigger
      GEE interpolation with no real precision gain.
    - Band: DSM (Digital Surface Model, metres above ellipsoid, EGM96)
"""

import sys
from pathlib import Path
import ee
import time

sys.path.append(str(Path(__file__).resolve().parents[1]))
from config.gee_config import GEE_PROJECT_ID, DEM_REGION

# ── GEE initialisation ────────────────────────────────────────────────────────
ee.Initialize(project=GEE_PROJECT_ID)

# ── Load and mosaic AW3D30 ────────────────────────────────────────────────────
dataset = ee.ImageCollection("JAXA/ALOS/AW3D30/V4_1").mosaic()

region = ee.Geometry.Rectangle(DEM_REGION)

elevation = dataset.select("DSM").clip(region)

# ── Export to Google Drive ────────────────────────────────────────────────────
task = ee.batch.Export.image.toDrive(
    image=elevation,
    description="Sihanoukville_ALOS_30m_DEM_V4_1",
    folder="GEE_Exports",
    fileNamePrefix="sihanoukville_alos30m_dem_1",
    region=region,
    scale=30,           # native resolution of AW3D30
    maxPixels=1e13,
    fileFormat="GeoTIFF",
)

task.start()
print("ALOS DEM export task submitted.")

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
