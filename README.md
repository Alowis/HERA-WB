# RegimeShift_codes

Pan-European hydrological regime shift analysis based on the HERA reanalysis (1951–2020).

This repository contains the analysis pipeline for:
- Deaggregation of nested catchment LISFLOOD outputs to residual (inter-catchment) values
- Validation against satellite products (GLEAM AET, ESA CCI soil moisture, GlobSnow SWE) and river gauges
- Trend detection (Sen's slope + Mann-Kendall) across ~2200 European catchments
- Water balance closure analysis
- Regime classification and attribution
- Publication-ready figure generation

## Repository Structure

```
RegimeShift_codes/
├── R/                          # Shared utility functions
├── scripts/
│   ├── 00_preprocessing/       # Data prep, deaggregation, aggregation
│   ├── 01_validation/          # Satellite & gauge comparison pipelines
│   ├── 02_trend_analysis/      # Monthly/yearly trend maps (all variables)
│   ├── 03_water_balance/       # Catchment & continental water balance
│   ├── 04_figures/             # Publication figures
│   ├── 05_regime_analysis/     # Regime classification & XAI
│   └── 06_export/              # CAMELS-style CSV & NetCDF export
├── python/                     # Python utilities (download, preprocessing)
├── tests/                      # Property-based tests
├── config/                     # User-editable path configuration
└── docs/                       # Variable descriptions & documentation
```

## Requirements

### R packages
```r
install.packages(c(
  "data.table", "sf", "terra", "exactextractr", "ncdf4",
  "ggplot2", "cowplot", "rnaturalearth", "rnaturalearthdata",
  "lubridate", "scales", "dplyr", "tidyr", "sp"
))
```

### Python packages
```
pip install pcraster netCDF4 numpy
```

### Data requirements (not included in repo)
- HERA LISFLOOD TSS output files (6-hourly, 1951–2020)
- Catchments GeoPackage (`catchments_analysis_final_v3.gpkg`)
- GLEAM v4.3a AET (daily/monthly NetCDF)
- GlobSnow v3.0 SWE (daily NetCDF)
- ESA CCI soil moisture
- Land use fraction NetCDFs (yearly, from LISFLOOD static maps)
- DEM and Köppen-Geiger climate classification rasters
- Observed river discharge (GRDC/national services)

## Quick Start

1. Edit `config/paths.R` to match your local data paths
2. Run preprocessing: `scripts/00_preprocessing/deaggregate_and_metrics.R`
3. Run validation pipelines in `scripts/01_validation/`
4. Generate trend maps: `scripts/02_trend_analysis/`
5. Check water balance: `scripts/03_water_balance/`
6. Produce figures: `scripts/04_figures/`

## Citation

Tilloy, A. et al. (2024). HERA: a high-resolution pan-European hydrological reanalysis (1951–2020).
Joint Research Centre, European Commission.

## License

[Add license]
