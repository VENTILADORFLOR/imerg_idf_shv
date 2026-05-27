# GPM IMERG IDF Correction — Sihanoukville

Reproducible code repository accompanying the paper:

> **Using Official IDF Benchmarks to Calibrate GPM IMERG Rainfall Intensities
> for Drainage Infrastructure Assessment in Data-Scarce Coastal Cambodia**
> Ying Li — School of Intelligent Construction, Wuchang University of Technology

**Study point:** 10.540139°N, 103.656389°E (Sihanoukville, Cambodia)  
**Satellite data:** NASA GPM IMERG V07 Final Run, half-hourly, 2000–2024  
**Official benchmark:** Sihanoukville IDF table, Department of Hydrology, 1999–2019

---

## Repository Structure

```
imerg_idf_sihanoukville/
│
├── config/
│   └── gee_config.py               ← GEE project ID, study point, DEM region
│
├── 01_data_extraction/             ← Step 1 — Extract raw data from GEE
│   ├── imerg_2000_2019.py          ← IMERG half-hourly series, 2000-06-01 to 2020-01-01
│   ├── imerg_2020_2024.py          ← IMERG half-hourly series, 2020-01-01 to 2025-01-01
│   └── get_alos_dem.py             ← ALOS AW3D30 DEM, 30 m, study region
│
├── 02_annual_maxima/               ← Step 2 — Rolling annual maxima
│   └── compute_annual_maxima.py    ← 5 durations × 2 periods; temporal comparison
│
├── 03_idf_fitting/                 ← Step 3 — Gumbel distribution fitting
│   └── gumbel_idf_2000_2019.py     ← IDF intensities for 5-, 10-, 25-year return periods
│
├── 04_watershed/                   ← Step 4 — Watershed delineation and figure
│   ├── delineate_watershed.py      ← Catchment boundary and flow path network (pysheds)
│   └── plot_catchment_figure.py    ← Publication-quality catchment map
│
├── 05_tables/                      ← Step 5 — Paper tables
│   └── generate_tables.py          ← Tables 1–4 from processed CSV outputs
│
├── data/
│   ├── raw/                        ← GEE exports + official IDF CSVs (not tracked by Git)
│   │   ├── rainfall_point_imerg_2000_2019.csv
│   │   ├── rainfall_point_imerg_2020_2024.csv
│   │   ├── sihanoukville_alos30m_dem.tif
│   │   ├── Design_Rainfall_Intensities_Selected.csv
│   │   └── IDF_TABLE.csv
│   ├── processed/                  ← Intermediate outputs (not tracked by Git)
│   ├── tables/                     ← Table inputs and outputs (not tracked by Git)
│   │   ├── Design_Rainfall_Intensities_Selected.csv  ← copy here before Step 5
│   │   └── IDF_TABLE.csv                             ← copy here before Step 5
│   ├── watershed/                  ← Shapefile outputs from Step 4
│   └── figures/                    ← Figure outputs from Step 4
│
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Workflow

### Prerequisites

```bash
pip install -r requirements.txt
earthengine authenticate          # one-time GEE authentication
```

Edit `config/gee_config.py` to set your own GEE project ID if needed.

---

### Step 1 — Extract raw data from Google Earth Engine

```bash
python 01_data_extraction/imerg_2000_2019.py
python 01_data_extraction/imerg_2020_2024.py
python 01_data_extraction/get_alos_dem.py
```

Tasks are submitted to GEE and run in the cloud (approximately 15–25 minutes each).
Download the exported files from Google Drive `GEE_Exports/` into `data/raw/`.

> **Note on IMERG V07 availability:**  
> IMERG V07 Final Run covers through September 2025, then production was suspended
> pending the V08 release (expected summer 2026). The 2020-2024 extraction uses
> `filterDate("2020-01-01", "2025-01-01")` — a right-open interval capturing
> five complete calendar years (2020–2024).

---

### Step 2 — Compute rolling annual maxima

```bash
python 02_annual_maxima/compute_annual_maxima.py
```

Outputs written to `data/processed/`:

| File | Description |
|------|-------------|
| `annual_maxima_imerg_2000_2019.csv` | Per-year maxima by duration, 2000–2019 |
| `annual_maxima_imerg_2020_2024.csv` | Per-year maxima by duration, 2020–2024 |
| `summary_2000_2019.csv` | Mean / max / min / std intensity per duration |
| `summary_2020_2024.csv` | Same for 2020–2024 period |
| `comparison_2000_2019_vs_2020_2024.csv` | Side-by-side with change % (paper Table 3) |
| `IMERG_temporal_stability_comparison.xlsx` | All of the above in one workbook |

---

### Step 3 — Gumbel distribution fitting

```bash
python 03_idf_fitting/gumbel_idf_2000_2019.py
```

Input: `data/processed/annual_maxima_imerg_2000_2019.csv`  
Output: `data/processed/imerg_idf_2000_2019_gumbel.csv`

Fits a two-parameter Gumbel distribution to each duration series and derives
design intensities for 5-, 10-, and 25-year return periods.

---

### Step 4 — Watershed delineation and catchment figure

```bash
python 04_watershed/delineate_watershed.py
python 04_watershed/plot_catchment_figure.py
```

`delineate_watershed.py` requires:

| File | Path |
|------|------|
| ALOS 30 m DEM | `data/raw/sihanoukville_alos30m_dem.tif` |
| Pour point | `data/drainage_node/drainage_node.shp` |

Outputs `data/watershed/catchment.shp` and `data/watershed/flow_paths.shp`.

> **Catchment area note:**  
> The Python (pysheds) delineation yields **3.599 km²**, compared to **3.43 km²**
> reported in the paper (derived from ArcGIS Pro). The difference (~5%) reflects
> two sources of discretisation error: (1) pysheds uses the native 30 m raster
> resolution, whereas the original ArcGIS analysis used a 12.5 m interpolated
> grid; (2) the D8 flow-direction algorithm implementation differs slightly
> between the two tools. The discrepancy does not affect the IDF correction
> methodology or conclusions.

> **DEM version note — AW3D30 V4.1 vs V3.2:**  
> The paper used V3.2; the extraction script has been updated to **V4.1**
> (released April 2024). V4.1 corrected anomalous elevation values in 19,051
> tiles globally using Copernicus DEM GLO-30 and ArcticDEM v4 as supplementary
> data. For the Sihanoukville study area (low-latitude coastal terrain), the
> elevation changes between V3.2 and V4.1 are negligible and do not affect the
> catchment delineation result. The GEE band name remains **`DSM`** (signed
> 16-bit integer, metres above ellipsoid, EGM96) — no script changes required.

`plot_catchment_figure.py` additionally requires:

| Layer | Path |
|-------|------|
| Site boundary | `data/study_area/study_area.shp` |
| Depression polygons | `data/depression/depression.shp` |
| Main road | `data/main_road/main_road.shp` |
| OSM road network | `data/roads/gis_osm_roads_free_1.shp` |

All shapefiles are reprojected to **EPSG:32648** (WGS 1984 UTM Zone 48N) at
load time, regardless of their source CRS.  
Output: `data/figures/catchment_figure.png`

---

### Step 5 — Generate paper tables

Copy the official IDF CSVs into `data/tables/`, then run:

```bash
python 05_tables/generate_tables.py
```

Inputs required:

| File | Description |
|------|-------------|
| `data/tables/Design_Rainfall_Intensities_Selected.csv` | Official benchmark intensities |
| `data/processed/imerg_idf_2000_2019_gumbel.csv` | IMERG Gumbel intensities (Step 3) |
| `data/processed/comparison_2000_2019_vs_2020_2024.csv` | Temporal comparison (Step 2) |

Outputs written to `data/tables/`:

| File | Paper table |
|------|-------------|
| `idf_comparison.csv` | Table 1 — Official vs IMERG intensities |
| `correction_factors.csv` | Table 2 — Duration correction factors αd |
| `temporal_stability.csv` | Table 3 — Temporal stability 2000–2019 vs 2020–2024 |
| `peak_discharge.csv` | Table 4 — Rational method peak discharge |
| `all_tables.xlsx` | All four tables in one workbook |

---

## Key Findings

| Duration | Correction factor αd |
|----------|--------------------|
| 30 min | 3.65 |
| 1 h | 2.90 |
| 2 h | 2.08 |
| 4 h | 1.69 |
| 24 h | 1.24 |

Raw IMERG underestimates short-duration design rainfall by **63–74%** relative
to the official IDF benchmark at 30–60 min durations. After correction, the
10-year 1-hour design intensity increases from 37.8 mm/h to 109.8 mm/h,
raising the estimated peak discharge from 19.8 m³/s to 57.6 m³/s.

---

## Data Sources

| Dataset | Source | Access |
|---------|--------|--------|
| GPM IMERG V07 Final Run | NASA / JAXA | Google Earth Engine `NASA/GPM_L3/IMERG_V07` |
| ALOS AW3D30 V4.1 DEM | JAXA | Google Earth Engine `JAXA/ALOS/AW3D30/V4_1` |
| OSM road network | Geofabrik Cambodia extract | [download.geofabrik.de](https://download.geofabrik.de) |
| Official IDF table | Department of Hydrology, Cambodia | Provided in `data/raw/` |

---

## Reproducibility Note

All IMERG and DEM data are static archives. Re-running the Step 1 extraction
scripts with the same parameters will produce byte-identical outputs.
The 2020–2024 IMERG series covers five complete calendar years; no partial-year
data is included, ensuring that annual maxima are comparable across years.
