# config/gee_config.py
# Central configuration for Google Earth Engine.
# All extraction scripts import from here — change project ID in one place only.

GEE_PROJECT_ID = "GEE_PROJECT_ID"

# Study point: Sihanoukville grid cell centroid (lon, lat)
STUDY_POINT = [103.656389, 10.540139]

# DEM export region bounding box [west, south, east, north]
DEM_REGION = [103.416998, 10.392559, 103.950503, 10.874039]
