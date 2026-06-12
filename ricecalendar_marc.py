def ricecalendar_marc(domain_path, basepath, referenceraster_path):
    """
    Download and preprocess the MARC (Map of Rice Cropping) rice calendar dataset
    (Laborte et al. 2017; https://db.cger.nies.go.jp/nies_data/10.17595/20230728.001/).

    Produces a NetCDF file with 12 layers:
        3 rice growing seasons (PaddyRice1, PaddyRice2, PaddyRice3)
        × 2 variables (planting DOY, growing_season_length)
        × 2 irrigation techniques (rf = rainfed, ir = irrigated)

    Note: MARC does not distinguish rainfed/irrigated — rf and ir layers are identical.
    """

    import os
    import glob
    import numpy as np
    import geopandas as gpd
    import xarray as xr
    import rioxarray  # noqa: F401  – registers .rio accessor on xarray objects
    from preproc_tools import download_url, unzip_all, ensure_xy_dims, makedirs, safe_clip

    # ------------------------------------------------------------------ #
    # Load reference raster and domain polygon
    # ------------------------------------------------------------------ #
    mask = gpd.read_file(domain_path)
    to_match = xr.open_dataset(referenceraster_path)
    to_match.rio.write_crs(4326, inplace=True)

    # ------------------------------------------------------------------ #
    # Download MARC rice calendar data
    # ------------------------------------------------------------------ #
    print("        *** DOWNLOADING MARC RICE CALENDAR DATA ***")
    url = 'https://db.cger.nies.go.jp/nies_data/10.17595/20230728.001/Rice_Calendar_Data_ver.1.0.zip'

    # Build target directory: rawdata/cropcalendar/marc_rice/
    cropcal_dir = makedirs(basepath, 'rawdata', 'cropcalendar')
    target_dir = os.path.join(cropcal_dir, 'marc_rice')
    os.makedirs(target_dir, exist_ok=True)

    download_path = os.path.join(target_dir, 'Rice_Calendar_Data_ver.1.0.zip')

    # Check for any already-unzipped content (NetCDF files present)
    existing_nc = glob.glob(os.path.join(target_dir, '**', '3_Transplanting_Harvest_Cropping.nc'), recursive=True)

    if existing_nc:
        nc_file = existing_nc[0]
        print(f"MARC rice calendar already unzipped, skipping: {nc_file}")
    else:
        if os.path.exists(download_path):
            print(f"ZIP already exists, skipping download: {download_path}")
        else:
            print("Downloading MARC rice calendar data")
            print('URL:', url)
            download_url(url, download_path=download_path)

        print("Unzipping MARC rice calendar data...")
        unzip_all(dir=target_dir)

        existing_nc = glob.glob(os.path.join(target_dir, '**', '3_Transplanting_Harvest_Cropping.nc'), recursive=True)
        if not existing_nc:
            raise FileNotFoundError(
                "Could not find '3_Transplanting_Harvest_Cropping.nc' after unzipping. "
                f"Check contents of: {target_dir}"
            )
        nc_file = existing_nc[0]

    # ------------------------------------------------------------------ #
    # Open source file
    # ------------------------------------------------------------------ #
    src = xr.open_dataset(nc_file)
    src = ensure_xy_dims(src)
    src.rio.write_crs(4326, inplace=True)

    # ------------------------------------------------------------------ #
    # Check that domain polygon lies fully within the raster extent
    # ------------------------------------------------------------------ #
    raster_bounds = src.rio.bounds()          # (left, bottom, right, top)
    poly_bounds = mask.to_crs(4326).total_bounds  # [minx, miny, maxx, maxy]

    fully_inside = (
        poly_bounds[0] >= raster_bounds[0] and   # poly west  ≥ raster west
        poly_bounds[1] >= raster_bounds[1] and   # poly south ≥ raster south
        poly_bounds[2] <= raster_bounds[2] and   # poly east  ≤ raster east
        poly_bounds[3] <= raster_bounds[3]        # poly north ≤ raster north
    )

    if fully_inside:
        print("Input polygon is fully inside the MARC raster extent.")
    else:
        print("WARNING: Input polygon extends (partly or fully) beyond the MARC raster extent!")
        print(f"  Polygon bounds [W, S, E, N]: {poly_bounds}")
        print(f"  Raster bounds  [W, S, E, N]: "
              f"[{raster_bounds[0]}, {raster_bounds[1]}, {raster_bounds[2]}, {raster_bounds[3]}]")

    # ------------------------------------------------------------------ #
    # Helper: iteratively fill NaN cells by copying the (plant, harvest)
    # pair from the valid 8-neighbor with the shortest growing season.
    # Both arrays are filled jointly so the pair always stays coherent.
    # ------------------------------------------------------------------ #
    def fill_nodata_joint(plant_arr, harvest_arr):
        """
        Fill NaN cells in plant_arr and harvest_arr jointly.
        For each NaN cell, examine up to 8 immediate neighbors that have BOTH
        valid planting and harvest values.  Among those valid neighbors, copy
        the (plant, harvest) pair from the one with the shortest growing season
        (harvest − plant, accounting for year-boundary crossing).
        Iterates until no further progress is possible.
        """
        offsets = [(-1, -1), (-1, 0), (-1, 1),
                   (0,  -1),          (0,  1),
                   (1,  -1), (1,  0), (1,  1)]

        plant   = plant_arr.copy().astype(float)
        harvest = harvest_arr.copy().astype(float)
        nrows, ncols = plant.shape
        prev_nan_count = -1

        while True:
            valid    = ~(np.isnan(plant) | np.isnan(harvest))
            nan_mask = ~valid
            nan_count = int(nan_mask.sum())
            if nan_count == 0 or nan_count == prev_nan_count:
                break   # fully filled, or no progress (isolated islands)

            # Season length for valid cells; invalid cells get inf so they
            # are never chosen as the "shortest-season" neighbor.
            raw_season = harvest - plant
            season = np.where(
                valid,
                np.where(raw_season <= 0, raw_season + 365, raw_season),
                np.inf
            )

            # Pad all three arrays to handle grid boundaries cleanly.
            pad_p = np.pad(plant,   1, mode='constant', constant_values=np.nan)
            pad_h = np.pad(harvest, 1, mode='constant', constant_values=np.nan)
            pad_s = np.pad(season,  1, mode='constant', constant_values=np.inf)

            # Build 8 shifted views: shape (8, nrows, ncols)
            nbr_p = np.stack([pad_p[1+dr:nrows+1+dr, 1+dc:ncols+1+dc] for dr, dc in offsets])
            nbr_h = np.stack([pad_h[1+dr:nrows+1+dr, 1+dc:ncols+1+dc] for dr, dc in offsets])
            nbr_s = np.stack([pad_s[1+dr:nrows+1+dr, 1+dc:ncols+1+dc] for dr, dc in offsets])

            # For each cell pick the neighbor index with the shortest season.
            best_idx     = np.argmin(nbr_s, axis=0)        # shape (nrows, ncols)
            best_idx_exp = best_idx[np.newaxis]             # shape (1, nrows, ncols)
            best_p = np.take_along_axis(nbr_p, best_idx_exp, axis=0)[0]
            best_h = np.take_along_axis(nbr_h, best_idx_exp, axis=0)[0]

            # Only fill NaN cells that have at least one valid neighbor.
            has_valid_nbr = np.min(nbr_s, axis=0) < np.inf
            fill_here = nan_mask & has_valid_nbr

            plant   = np.where(fill_here, best_p, plant)
            harvest = np.where(fill_here, best_h, harvest)
            prev_nan_count = nan_count

        return plant, harvest

    # ------------------------------------------------------------------ #
    # Define crop cycles
    # Cycle 1 → PaddyRice1, cycle 2 → PaddyRice2, cycle 3 → PaddyRice3
    # Variable name 'Transplating' preserves the original (misspelled) name
    # in the MARC NetCDF file.
    # ------------------------------------------------------------------ #
    crop_cycles = [
        ('PaddyRice1', 'Transplanting_1_Cropping', 'Harvest_1_Cropping'),
        ('PaddyRice2', 'Transplanting_2_Cropping', 'Harvest_2_Cropping'),
        ('PaddyRice3', 'Transplanting_3_Cropping', 'Harvest_3_Cropping'),
    ]

    file_to_mosaic = []

    for crop_name, transplant_var, harvest_var in crop_cycles:

        print(f"    Processing {crop_name} ({transplant_var} / {harvest_var})")

        # Extract DataArrays; xarray auto-masks fill values as NaN
        plant_da = src[transplant_var].astype(float)
        harvest_da = src[harvest_var].astype(float)

        # Treat zero or negative DOY as nodata (safety guard for non-CF fill values)
        plant_da = plant_da.where(plant_da > 0)
        harvest_da = harvest_da.where(harvest_da > 0)

        # ---- Fill nodata BEFORE resampling (joint fill, shortest-season neighbor) ----
        plant_filled_arr, harvest_filled_arr = fill_nodata_joint(plant_da.values, harvest_da.values)

        plant_da_filled   = plant_da.copy(data=plant_filled_arr)
        harvest_da_filled = harvest_da.copy(data=harvest_filled_arr)

        # Propagate CRS to the filled DataArrays
        plant_da_filled.rio.write_crs(4326, inplace=True)
        harvest_da_filled.rio.write_crs(4326, inplace=True)

        # ---- Compute growing season length (cyclic DOY) ----
        season_days = harvest_da_filled - plant_da_filled
        # When harvest DOY < planting DOY the season crosses the year boundary
        season_days = xr.where(season_days <= 0, season_days + 365, season_days)
        season_days.rio.write_crs(4326, inplace=True)

        # ---- Reproject to reference grid ----
        plant_proj = plant_da_filled.rio.reproject_match(to_match)
        season_proj = season_days.rio.reproject_match(to_match)

        # ---- Clip to domain polygon ----
        plant_clip = safe_clip(plant_proj, mask)
        season_clip = safe_clip(season_proj, mask)

        # Replace DOY 60 (Feb 29) with 61 (Mar 1) to avoid leap-year artefacts
        plant_clip = plant_clip.where(plant_clip != 60, 61)

        # ---- Create rf and ir layers (identical — MARC has no irrigation split) ----
        for technique in ['rf', 'ir']:
            layername = f'{crop_name}_{technique}'

            plant_layer = xr.DataArray(
                np.round(plant_clip.values).astype('float32'),
                dims=plant_clip.dims,
                coords=plant_clip.coords,
                name=f'{layername}_planting',
                attrs={'units': 'day of year',
                       'long_name': f'Planting day {technique}'},
            )

            # Store as float32 with integer values (range ~25–200) so that the
            # values are directly readable without timedelta decoding.
            season_layer = xr.DataArray(
                np.round(season_clip.values).astype('float32'),
                dims=season_clip.dims,
                coords=season_clip.coords,
                name=f'{layername}_growing_season_length',
                attrs={'units': 'days',
                       'long_name': f'Growing season length {technique}'},
            )

            file_to_mosaic.append(plant_layer)
            file_to_mosaic.append(season_layer)

    # ------------------------------------------------------------------ #
    # Merge all layers and save
    # ------------------------------------------------------------------ #
    src_mosaic = xr.merge(file_to_mosaic)
    src_mosaic = src_mosaic.rio.write_crs(4326)   # adds spatial_ref coordinate

    target_out_dir = makedirs(basepath, 'processed', 'cropcalendar_marc')
    targetfile = os.path.join(target_out_dir, 'cropcalendar.nc')
    src_mosaic.to_netcdf(targetfile)

    print(f"MARC rice calendar saved to: {targetfile}")
    print(f"Output variables ({len(src_mosaic.data_vars) - 1} + spatial_ref):")
    for var in src_mosaic.data_vars:
        if var != 'spatial_ref':
            print(f"  {var}")

domain_path = '/Users/ritterj1/PythonProjects/aquacropgrid-preproc/inputdata/mekong/basin_outline/mekong_jrc_outline.geojson'
basepath = '/Users/ritterj1/PythonProjects/aquacropgrid-preproc'
referenceraster_path = '/Users/ritterj1/PythonProjects/aquacropgrid-preproc/template_grid.nc'
ricecalendar_marc(domain_path, basepath, referenceraster_path)