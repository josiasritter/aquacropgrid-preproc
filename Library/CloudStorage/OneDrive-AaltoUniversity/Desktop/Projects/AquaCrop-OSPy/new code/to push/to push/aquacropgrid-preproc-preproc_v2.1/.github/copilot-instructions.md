# Copilot Instructions for AquaCropGrid-Preproc

## Project Overview
- This project preprocesses global datasets (climate, soil, crop areas, crop calendars) for use with AquaCrop-OSPy in a gridded format.
- Main workflow: Download, preprocess, and format input data into NetCDF for model input.
- Key files:
  - `preproc_main.py`: Entry point; configure input arguments here before running.
  - `preproc_tools.py`: Core preprocessing, download, and utility functions.
  - `cropcalendar_module.py`, `climate.py`, `soil.py`, `crop_areas.py`: Specialized modules for each data type.
  - `environment.yml`: Conda environment definition (Python 3.10, see dependencies).

## Developer Workflow
- **Environment setup:**
  - Run `conda env create -f environment.yml` and `conda activate aquacropgrid-preproc`.
- **Running the pipeline:**
  - Edit input arguments at the top of `preproc_main.py` (paths, years, etc.).
  - Run: `python preproc_main.py`
- **Data sources:**
  - Some datasets require manual download or API tokens (see README for Copernicus CDS setup).
  - Downloaded and processed data is organized under `rawdata/`, `inputdata/`, and `processed/`.

## Project Conventions & Patterns
- All geospatial processing uses `xarray`, `rioxarray`, `geopandas`, and `rasterio`.
- Helper functions for downloading/unzipping are in `preproc_tools.py`.
- NetCDF output follows CF-1.8 conventions; see `basegrid()` in `preproc_tools.py` for template.
- Masking, reprojection, and clipping are handled via `geopandas` and `rioxarray`.
- Crop calendar mapping uses a dictionary to translate GGCMI crop IDs to AquaCrop names.
- Avoid re-downloading or re-unzipping if files already exist (see checks in modules).

## Integration & Extensibility
- To add new data sources, create a new module and follow the structure in `preproc_tools.py`.
- Use `makedirs()` and `ensure_xy_dims()` helpers for consistent directory and data handling.
- All scripts expect relative paths from the project root.

## Testing & Debugging
- No formal test suite; validate outputs by inspecting generated NetCDF files in `processed/`.
- Use `pdb.set_trace()` for debugging (already imported in `preproc_tools.py`).
- For troubleshooting, check for missing dependencies or incorrect file paths.

## Example: Adding a New Preprocessing Step
1. Add a function to `preproc_tools.py` or a new module.
2. Import and call it from `preproc_main.py`.
3. Use existing helpers for file I/O and geospatial operations.

---
For more details, see `README.md` and code comments in each module.
