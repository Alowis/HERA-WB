# HERA-WB

**Hydrological European ReAnalysis – Water Balance** (1951–2020)

This repository contains the analysis pipeline for evaluating and exporting water balance components from the HERA LISFLOOD reanalysis at the catchment scale across Europe.

The pipeline covers:
- Deaggregation of nested catchment LISFLOOD outputs to residual (inter-catchment) values
- Validation against benchmark products (GLEAM AET, ESA CCI soil moisture, GlobSnow SWE) and river discharge gauges
- Water balance closure analysis
- CAMELS-style data export (CSV and Parquet)
- Publication-ready figure generation

## Repository Structure

```
HERA-WB/
├── config/
│   └── paths.R                       # User-editable path configuration
├── docs/
│   └── variable_descriptions.csv     # Metadata for all HERA variables
├── R/
│   ├── build_workspace.R             # Build shared .rds workspace (run once)
│   ├── load_workspace.R              # Loader sourced by all downstream scripts
│   └── functions_regime.R            # Utility functions (regime, NetCDF I/O, Budyko)
├── python/
│   ├── download_gleam_daily.py       # Download GLEAM v4.3a daily AET via SFTP
│   ├── download_globsnow.py          # Download GlobSnow SWE
│   └── LF_tssPreProcess.py          # LISFLOOD TSS pre-processing (PCRaster)
├── scripts/
│   ├── 00_preprocessing/
│   │   ├── 0.Aggregate_landuse_catchments.R
│   │   ├── 0.Compute_catchment_attributes.R
│   │   ├── 1.Nested_Catchments.R
│   │   ├── 2.Deaggregate_and_metrics.R
│   │   └── 2.Preprocess_tss_aggregates.R
│   ├── 01_evaluation/
│   │   ├── aet/                      # GLEAM AET validation
│   │   ├── dis/                      # River discharge validation
│   │   ├── ssm/                      # ESA CCI soil moisture validation
│   │   └── swe/                      # GlobSnow SWE validation
│   ├── 02_figures/
│   │   ├── Fig01_nesting_map.R
│   │   ├── Fig02_combined_AET_SM.R
│   │   ├── Fig03_combined_SWE_discharge.R
│   │   ├── Fig04_evaluation_summary.R
│   │   ├── Fig05_panels.R
│   │   └── Supplement/               # Supplementary figures
│   └── 03_export/
│       ├── export_camels_style.R     # CAMELS-style CSV per catchment
│       └── export_landuse_camels_style.R
│       └── export_parquet.R

└── HERA-WB.Rproj                     # RStudio project file
```

## Requirements

### R packages

```r
install.packages(c(
  "data.table", "sf", "terra", "exactextractr", "ncdf4",
  "ggplot2", "cowplot", "rnaturalearth", "rnaturalearthdata",
  "lubridate", "scales", "dplyr", "tidyr", "sp",
  "raster", "viridis", "RtsEva", "moments", "ineq"
))
```

### Python packages

```
pip install paramiko pcraster netCDF4 numpy pandas
```

### Data requirements (not included in repo)

- HERA LISFLOOD TSS output files (6-hourly, 1951–2020)
- Catchments GeoPackage (`catchments_analysis_final_v3.gpkg`)
- GLEAM v4.3a AET (daily/monthly NetCDF)
- GlobSnow v3.0 SWE (daily NetCDF)
- ESA CCI soil moisture
- Land use fraction NetCDFs (yearly, from LISFLOOD static maps)
- DEM and Köppen-Geiger climate classification rasters
- Observed river discharge (GRDC / national services)

## Quick Start

1. **Configure paths** — Edit `config/paths.R` to point to your local data directories.

2. **Download satellite products** (if needed):
   ```bash
   python python/download_gleam_daily.py
   python python/download_globsnow.py
   ```

3. **Pre-process LISFLOOD outputs**:
   ```bash
   python python/LF_tssPreProcess.py
   ```

4. **Run preprocessing scripts** in order:
   - `scripts/00_preprocessing/0.Aggregate_landuse_catchments.R`
   - `scripts/00_preprocessing/0.Compute_catchment_attributes.R`
   - `scripts/00_preprocessing/1.Nested_Catchments.R`
   - `scripts/00_preprocessing/2.Deaggregate_and_metrics.R`
   - `scripts/00_preprocessing/2.Preprocess_tss_aggregates.R`

5. **Build the shared workspace** (run once after preprocessing):
   ```r
   source("R/build_workspace.R")
   ```
   This produces `output/workspace.rds`, which all downstream scripts load via `source("R/load_workspace.R")`.

6. **Run evaluation pipelines** — Scripts in `scripts/01_evaluation/` validate HERA outputs against satellite and gauge data (AET, soil moisture, SWE, discharge).

7. **Generate figures** — Scripts in `scripts/02_figures/` produce publication figures.

8. **Export data** — Scripts in `scripts/03_export/` export catchment time series and attributes in CAMELS-style CSV and NetCDF formats.

## Variables

The full list of HERA variables with units and descriptions is in [`docs/variable_descriptions.csv`](docs/variable_descriptions.csv). Key variables include:

| Variable | Description | Units |
|----------|-------------|-------|
| rainfall | Liquid precipitation | mm/6h |
| snowfall | Solid precipitation | mm/6h |
| ActEvapo | Actual evapotranspiration | mm/6h |
| theta1/2/3 | Soil moisture (upper/middle/deep) | m3/m3 |
| snowSWE | Snow water equivalent | mm |
| disWin | Channel discharge | m³/s |
| qlz | Baseflow (lower groundwater zone) | mm/6h |


## License

