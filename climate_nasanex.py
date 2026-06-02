"""
This script downloads and preprocesses the following daily climate projection data
from NASA NEX-GDDP-CMIP6:
    - Minimum temperature [°C]           (tasmin)
    - Maximum temperature [°C]           (tasmax)
    - Precipitation [mm/day]             (pr)
    - Reference evapotranspiration       computed from tasmin, tasmax, hurs,
      [mm/day, FAO-56 Penman-Monteith]   rsds, and sfcWind via the FAO-56 PM eq.

Data source: https://www.nccs.nasa.gov/data-collections/nex-gddp-cmip6/
Spatial subsets are fetched via the NCCS THREDDS NetCDF Subset Service (NCSS);
no authentication or AWS CLI is required.

How historical / SSP splitting works
--------------------------------------
NEX-GDDP-CMIP6 stores 1950-2014 under the 'historical' scenario and 2015-2100
under the chosen SSP. This script automatically fetches the correct scenario for
each year: years <= 2014 always use 'historical', years >= 2015 use the scenario
argument.  If scenario='historical' is passed explicitly, all years must be <= 2014.

Unit conversions applied
------------------------
    pr       : kg m⁻² s⁻¹  →  mm day⁻¹  (× 86400)
    tasmin   : K             →  °C        (− 273.15)
    tasmax   : K             →  °C        (− 273.15)
    rsds     : W m⁻²        →  MJ m⁻² day⁻¹  (× 0.0864)
    sfcWind  : m s⁻¹ at 10 m  →  m s⁻¹ at 2 m  (FAO-56 Eq. 47)

See https://ds.nccs.nasa.gov/thredds/catalog/AMES/NEX/GDDP-CMIP6/catalog.html
for the full list of available models and their ensemble members.
"""

import os
import time
import xml.etree.ElementTree as ET

import numpy as np
import requests
import scipy.ndimage as ndi
import xarray as xr
from rasterio.warp import Resampling

from preproc_tools import ensure_xy_dims, makedirs

# Cache of (model, scenario, ensemble, variable) -> (grid_label, version_suffix)
# populated lazily by _discover_file_pattern().
_pattern_cache: dict = {}


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

def _discover_file_pattern(model, scenario, ensemble, variable):
    """
    Query the THREDDS catalog XML to discover the grid label and optional
    version suffix used in the actual NetCDF filenames.

    Different models use different grid labels (e.g. 'gn' vs 'gr1') and some
    files carry a version suffix (e.g. '_v2.0').  This function fetches the
    folder listing once and caches the result so subsequent calls are free.

    Returns
    -------
    (grid, version_suffix) : (str, str)
        e.g. ('gr1', '_v2.0') for GFDL-CM4, or ('gn', '') for ACCESS-CM2.
    """
    key = (model, scenario, ensemble, variable)
    if key in _pattern_cache:
        return _pattern_cache[key]

    catalog_url = (
        "https://ds.nccs.nasa.gov/thredds/catalog/AMES/NEX/GDDP-CMIP6"
        f"/{model}/{scenario}/{ensemble}/{variable}/catalog.xml"
    )
    prefix = f"{variable}_day_{model}_{scenario}_{ensemble}_"

    try:
        resp = requests.get(catalog_url, timeout=30)
        resp.raise_for_status()
        root = ET.fromstring(resp.content)
        # Iterate all XML elements regardless of namespace
        for elem in root.iter():
            name = elem.get('name', '')
            if name.startswith(prefix) and name.endswith('.nc'):
                # Filename: {var}_day_{model}_{scenario}_{ensemble}_{grid}_{year}[_{version}].nc
                tail = name[len(prefix):-3]   # e.g. 'gr1_2015_v2.0'
                parts = tail.split('_')
                grid = parts[0]              # e.g. 'gr1'
                # parts[1] is the year; anything after is the version
                version_suffix = ('_' + '_'.join(parts[2:])) if len(parts) > 2 else ''
                _pattern_cache[key] = (grid, version_suffix)
                print(f"        Discovered filename pattern for {model}/{variable}: "
                      f"grid='{grid}', version='{version_suffix}'")
                return grid, version_suffix
    except Exception as exc:
        print(f"        Warning: THREDDS catalog query failed ({exc}). "
              f"Defaulting to grid='gn', no version suffix.")

    _pattern_cache[key] = ('gn', '')
    return 'gn', ''


def _scenario_for_year(year, ssp):
    """Return 'historical' for years up to 2014, otherwise the given SSP."""
    return 'historical' if year <= 2014 else ssp


def _build_ncss_url(model, scenario, ensemble, variable, year, north, west, south, east,
                    grid='gn', version_suffix=''):
    """Build a THREDDS NCSS spatial-subset download URL."""
    base = "https://ds.nccs.nasa.gov/thredds/ncss/grid/AMES/NEX/GDDP-CMIP6"
    fname = f"{variable}_day_{model}_{scenario}_{ensemble}_{grid}_{year}{version_suffix}.nc"
    path = f"{base}/{model}/{scenario}/{ensemble}/{variable}/{fname}"
    params = (
        f"?var={variable}"
        f"&north={north}&west={west}&east={east}&south={south}"
        "&horizStride=1"
        f"&time_start={year}-01-01T12:00:00Z"
        f"&time_end={year}-12-31T12:00:00Z"
        "&accept=netcdf4"
    )
    return path + params


def _download_year(model, ssp, ensemble, variable, year,
                   north, west, south, east, target_dir,
                   max_retries=5):
    """
    Download one spatially-subsetted yearly file from THREDDS NCSS.

    Returns the path to the local file (already downloaded or just saved).
    Retries up to *max_retries* times on transient network errors.
    """
    scenario = _scenario_for_year(year, ssp)
    local_fname = f"{variable}_{model}_{scenario}_{year}.nc"
    local_path = os.path.join(target_dir, local_fname)

    if os.path.exists(local_path):
        return local_path

    # Discover the correct grid label and version suffix for this model
    grid, version_suffix = _discover_file_pattern(model, scenario, ensemble, variable)

    url = _build_ncss_url(model, scenario, ensemble, variable, year,
                          north, west, south, east, grid, version_suffix)
    print(f"        *** DOWNLOADING NASA NEX-GDDP-CMIP6: "
          f"{model} {scenario} {variable} {year} ***")

    for attempt in range(1, max_retries + 1):
        try:
            response = requests.get(url, timeout=180, stream=True)
            response.raise_for_status()
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=1 << 20):  # 1 MB
                    f.write(chunk)
            return local_path
        except Exception as exc:
            if attempt < max_retries:
                wait = 15 * attempt
                print(f"        Attempt {attempt} failed ({exc}). "
                      f"Retrying in {wait} s…")
                time.sleep(wait)
            else:
                raise RuntimeError(
                    f"Could not download {url} after {max_retries} attempts."
                ) from exc


# ---------------------------------------------------------------------------
# FAO-56 Penman-Monteith ET₀  (vectorised over xarray)
# ---------------------------------------------------------------------------

def _calc_et0_xr(tasmin, tasmax, hurs, rsds, wind10m, elev=0.0):
    """
    Compute daily FAO-56 Penman-Monteith reference ET₀ from xarray DataArrays.

    All spatial DataArrays must share the same (time, y, x) dimensions; y must
    represent latitude in degrees (EPSG:4326).

    Parameters
    ----------
    tasmin, tasmax : xr.DataArray  – near-surface temperature min/max [°C]
    hurs           : xr.DataArray  – near-surface relative humidity [%]
    rsds           : xr.DataArray  – surface downwelling shortwave [MJ m⁻² day⁻¹]
    wind10m        : xr.DataArray  – wind speed at 10 m height [m s⁻¹]
    elev           : float         – mean domain elevation [m a.s.l.]; default 0

    Returns
    -------
    et0 : xr.DataArray [mm day⁻¹], non-negative, named 'ReferenceET'
    """
    # --- Physical constants --------------------------------------------------
    Gsc  = 0.0820    # solar constant              [MJ m⁻² min⁻¹]
    sbc  = 4.903e-9  # Stefan-Boltzmann constant   [MJ m⁻² day⁻¹ K⁻⁴]
    alb  = 0.23      # grass albedo (FAO-56)
    Cn   = 900.0     # short-crop numerator constant
    Cd   = 0.34      # short-crop denominator constant
    G    = 0.0       # soil heat flux ≈ 0 for daily time step

    # --- Atmospheric pressure and psychrometric constant ---------------------
    AtmP = 101.3 * ((293.0 - 0.0065 * elev) / 293.0) ** 5.26  # [kPa]
    psy  = 0.000665 * AtmP                                       # [kPa °C⁻¹]

    # --- Wind speed: 10 m → 2 m  (FAO-56 Eq. 47) ----------------------------
    u2 = wind10m * (4.87 / np.log(67.8 * 10.0 - 5.42))

    # --- Temperature and humidity --------------------------------------------
    tmean = (tasmax + tasmin) / 2.0
    rh    = hurs.clip(0.0, 100.0)

    # Saturation vapour pressure [kPa]
    e0max = 0.6108 * np.exp(17.27 * tasmax / (tasmax + 237.3))
    e0min = 0.6108 * np.exp(17.27 * tasmin / (tasmin + 237.3))
    es    = (e0max + e0min) / 2.0

    # Actual vapour pressure [kPa]
    ea = (rh / 100.0) * es

    # Slope of saturation vapour pressure curve [kPa °C⁻¹]
    Delta = (
        4098.0 * 0.6108 * np.exp(17.27 * tmean / (tmean + 237.3))
    ) / (tmean + 237.3) ** 2

    # --- Radiation -----------------------------------------------------------
    # Julian day of year (1-D along 'time')
    J = tasmin.time.dt.dayofyear.astype(float)

    # Latitude in radians (1-D along 'y')
    lat_rad = np.deg2rad(tasmin.y)

    # Inverse relative Earth-Sun distance and solar declination
    dr       = 1.0 + 0.033 * np.cos(2.0 * np.pi / 365.0 * J)
    sol_decl = 0.409 * np.sin(2.0 * np.pi / 365.0 * J - 1.39)

    # Sunset hour angle [rad]:  (-tan φ · tan δ) must be clamped to [-1, 1]
    # xarray broadcasts (y,) × (time,)  →  (time, y) automatically
    arccos_arg = (-np.tan(lat_rad) * np.tan(sol_decl)).clip(-1.0, 1.0)
    ws = np.arccos(arccos_arg)  # shape (time, y)

    # Extraterrestrial radiation [MJ m⁻² day⁻¹]  shape (time, y)
    Ra = (
        (24.0 * 60.0 / np.pi) * Gsc * dr * (
            ws * np.sin(lat_rad) * np.sin(sol_decl)
            + np.cos(lat_rad) * np.cos(sol_decl) * np.sin(ws)
        )
    )

    # Net shortwave radiation [MJ m⁻² day⁻¹]  shape (time, y, x)
    Rns = (1.0 - alb) * rsds

    # Clear-sky radiation [MJ m⁻² day⁻¹];  avoid division by zero (polar night)
    Rs0 = (0.75 + 2.0e-5 * elev) * Ra
    Rs0 = Rs0.where(Rs0 > 0.0, other=np.nan)

    # Cloudiness fraction, clamped to [0.05, 1.0]
    fcd = (1.35 * rsds / Rs0 - 0.35).clip(0.05, 1.0)
    fcd = fcd.where(Ra > 0.0, other=0.05)  # set to minimum during polar night

    # Net longwave radiation [MJ m⁻² day⁻¹]
    Rnl = (
        sbc * fcd
        * (0.34 - 0.14 * np.sqrt(ea.clip(0.0)))
        * (((tasmax + 273.16) ** 4 + (tasmin + 273.16) ** 4) / 2.0)
    )

    # Net radiation [MJ m⁻² day⁻¹]
    Rn = Rns - Rnl

    # --- FAO-56 PM ET₀ [mm day⁻¹] -------------------------------------------
    et0 = (
        0.408 * Delta * (Rn - G)
        + psy * (Cn / (tmean + 273.0)) * u2 * (es - ea)
    ) / (Delta + psy * (1.0 + Cd * u2))

    # Set ET₀ to 0 during polar night and floor at 0
    et0 = et0.where(Ra > 0.0, other=0.0).clip(0.0)

    et0.name = 'ReferenceET'
    et0.attrs.update({
        'long_name': 'FAO-56 Penman-Monteith reference evapotranspiration',
        'units': 'mm/day',
        'standard_name': 'water_potential_evapotranspiration',
    })
    return et0


# ---------------------------------------------------------------------------
# Preprocessing: reproject → gap-fill → mask → save
# ---------------------------------------------------------------------------

def _preproc_and_save(src, variable, yearlist, basepath, to_match, model, scenario, ensemble):
    """
    Reproject *src* to the project grid, gap-fill NaN, apply domain mask,
    and write the result to processed/.

    Parameters
    ----------
    src      : xr.Dataset with a single data variable named *variable*
    variable : str  AquaCrop variable name (e.g. 'MaxTemp')
    yearlist : list[int]
    basepath : str
    to_match : xr.Dataset  template raster from basegrid()
    model    : str  CMIP6 model identifier (stored as file attribute)
    scenario : str  SSP scenario identifier (stored as file attribute)
    ensemble : str  ensemble member identifier (stored as file attribute)
    """
    print(f"        *** PREPROCESSING NASA NEX-GDDP-CMIP6: {variable} ***")

    src = ensure_xy_dims(src)

    # Convert cftime/object time to proper datetime64[ns] and normalise to midnight
    try:
        src = src.assign_coords(
            time=src.indexes['time'].to_datetimeindex().normalize()
        )
    except (AttributeError, TypeError):
        pass  # already datetime64[ns]

    src.rio.write_crs(4326, inplace=True)

    # Resample to project grid
    src_reproj = src.rio.reproject_match(to_match, resampling=Resampling.nearest)

    # Spatial gap-fill: nearest-neighbour interpolation for ocean / data voids
    data_vars = list(src_reproj.data_vars)
    for vname in data_vars:
        data3d = src_reproj[vname].to_numpy()
        # Use the first time step to identify consistently missing cells
        nodata_mask = np.isnan(data3d[0])
        if nodata_mask.any():
            _, nearest_idx = ndi.distance_transform_edt(
                nodata_mask, return_indices=True
            )
            src_reproj[vname].data = data3d[:, nearest_idx[0], nearest_idx[1]]

    # Apply domain mask (cells outside polygon → NaN)
    src_masked = src_reproj.where(to_match['Band1'] == 1)

    # Build output path
    target_dir = makedirs(basepath, 'processed', '')
    out_name = f"{variable}{yearlist[0]}{yearlist[-1]}.nc"
    targetfile = os.path.join(target_dir, out_name)

    # Drop singleton dims / auxiliary coords; re-attach CRS
    src_masked = src_masked.squeeze(drop=True)
    src_masked = src_masked.drop_vars('spatial_ref', errors='ignore')
    src_masked = src_masked.rio.write_crs(4326)

    # Cast all data variables to float32 (matches AgERA5; ET0 computation produces float64)
    for vname in data_vars:
        if src_masked[vname].dtype != np.float32:
            src_masked[vname] = src_masked[vname].astype(np.float32)

    # Standardised variable attributes aligned with AgERA5 conventions
    _var_attrs = {
        'Precipitation': {
            'units': 'mm d-1', 'long_name': 'Total precipitation',
            'temporal_aggregation': 'Sum 00-00LT',
        },
        'ReferenceET': {
            'units': 'mm d-1',
            'long_name': 'Penman-Monteith reference evapotranspiration according to the FAO56 approach',
            'temporal_aggregation': 'Sum 00-00LT',
        },
        'MinTemp': {'units': 'degC', 'long_name': 'Minimum daily air temperature'},
        'MaxTemp': {'units': 'degC', 'long_name': 'Maximum daily air temperature'},
    }
    for vname in data_vars:
        src_masked[vname].encoding.pop('grid_mapping', None)
        attrs = dict(src_masked[vname].attrs)
        attrs.update(_var_attrs.get(variable, {}))
        attrs['grid_mapping'] = 'spatial_ref'
        src_masked[vname].attrs = attrs

    # Global CF attributes + NASA NEX provenance
    import datetime
    src_masked.attrs = {
        'Conventions' : 'CF-1.7',
        'title'       : 'NASA NEX-GDDP-CMIP6 downscaled daily climate projections.',
        'source'      : f'NASA NEX-GDDP-CMIP6 – model: {model}, scenario: {scenario}, ensemble: {ensemble}',
        'cmip6_model'    : model,
        'cmip6_scenario' : scenario,
        'cmip6_ensemble' : ensemble,
        'history'     : f'Preprocessed on {datetime.date.today().isoformat()} by aquacropgrid-preproc.',
        'references'  : 'https://www.nasa.gov/nex/gddp',
    }

    src_masked.to_netcdf(
        targetfile, mode='w',
        encoding={data_vars[0]: {'zlib': True, 'complevel': 4}},
    )


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def climate_nasanex(
    basepath,
    start_year,
    end_year,
    to_match,
    model='GFDL-CM4',
    scenario='ssp245',
    ensemble='r1i1p1f1',
    variables=None,
    elev=0.0,
):
    """
    Download and preprocess NASA NEX-GDDP-CMIP6 daily climate projections for
    a given area of interest and time period.

    Parameters
    ----------
    basepath   : str
        Working directory (must contain template_grid.nc).
    start_year : int
        First year to process.  Historical range: 1950–2014;
        SSP range: 2015–2100.  Years across both ranges are handled
        automatically.
    end_year   : int
        Last year to process (inclusive).
    to_match   : xr.Dataset
        Template raster returned by preproc_tools.basegrid(); defines the
        output grid and domain mask.
    model      : str
        CMIP6 model name.  Default 'GFDL-CM4'.
        Full list: https://ds.nccs.nasa.gov/thredds/catalog/AMES/NEX/GDDP-CMIP6/catalog.html
        Commonly used models include: ACCESS-CM2, BCC-CSM2-MR, CESM2,
        CMCC-CM2-SR5, GFDL-CM4, GFDL-ESM4, GISS-E2-1-G, MRI-ESM2-0.
    scenario   : str
        Greenhouse gas scenario for years >= 2015.
        Options: 'ssp126', 'ssp245', 'ssp370', 'ssp585'.
        Years <= 2014 always use the 'historical' scenario automatically.
        Pass 'historical' explicitly if only downloading historical data.
    ensemble   : str
        Ensemble member identifier.  Default 'r1i1p1f1'.
        Check the THREDDS catalog for the correct member for your model.
    variables  : list[str] or None
        Subset of ['MinTemp', 'MaxTemp', 'Precipitation', 'ReferenceET'].
        Defaults to all four when None.
    elev       : float
        Representative elevation of the domain [m a.s.l.], used in the ET₀
        atmospheric pressure term.  Default 0 m (sea level).

    Output files
    ------------
    Preprocessed NetCDF files are written to <basepath>/processed/ with names:
        {variable}_{model}_{scenario}_{start_year}{end_year}.nc
    Raw downloaded files are kept in <basepath>/rawdata/climate_nasanex/.
    """
    if variables is None:
        variables = ['MinTemp', 'MaxTemp', 'Precipitation', 'ReferenceET']

    yearlist   = list(range(start_year, end_year + 1))
    target_dir = makedirs(basepath, 'rawdata', 'climate_nasanex')

    # Bounding box from template grid (EPSG:4326, -180/180 convention)
    to_match.rio.write_crs(4326, inplace=True)
    xmin, ymin, xmax, ymax = [round(b, 3) for b in to_match.rio.bounds()]
    west, east, south, north = xmin, xmax, ymin, ymax

    # Map AquaCrop names → NEX-GDDP-CMIP6 variable names
    _direct = {
        'MinTemp'     : 'tasmin',
        'MaxTemp'     : 'tasmax',
        'Precipitation': 'pr',
    }
    # Variables required to compute ET₀ via FAO-56 PM
    _et0_inputs = ['tasmin', 'tasmax', 'hurs', 'rsds', 'sfcWind']

    # Determine which raw NEX variables need downloading
    nex_vars_needed = set()
    for v in variables:
        if v in _direct:
            nex_vars_needed.add(_direct[v])
        elif v == 'ReferenceET':
            nex_vars_needed.update(_et0_inputs)

    # ------------------------------------------------------------------
    # Step 1 – Download: one spatially-subsetted file per variable/year
    # ------------------------------------------------------------------
    for nex_var in sorted(nex_vars_needed):
        for year in yearlist:
            _download_year(
                model, scenario, ensemble,
                nex_var, year,
                north, west, south, east,
                target_dir,
            )

    # ------------------------------------------------------------------
    # Step 2 – Helper: load and concatenate yearly files for one NEX var
    # ------------------------------------------------------------------
    def _load_years(nex_var):
        files = []
        for year in yearlist:
            scen = _scenario_for_year(year, scenario)
            files.append(
                os.path.join(target_dir, f"{nex_var}_{model}_{scen}_{year}.nc")
            )
        datasets = [xr.open_dataset(f) for f in files]
        ds = xr.concat(datasets, dim='time').sortby('time')
        # Normalise any 0-360 longitudes to -180/180
        ds = ensure_xy_dims(ds)
        if 'x' in ds.coords and float(ds.x.max()) > 180.0:
            ds = ds.assign_coords(x=((ds.x + 180.0) % 360.0 - 180.0)).sortby('x')
        return ds

    # ------------------------------------------------------------------
    # Step 3 – Process each requested output variable
    # ------------------------------------------------------------------
    for variable in variables:

        if variable in ('MinTemp', 'MaxTemp'):
            nex_var = _direct[variable]
            ds = _load_years(nex_var)
            da = ds[nex_var] - 273.15          # K → °C
            da.attrs.update({'units': 'degC'})
            out = da.to_dataset(name=variable)

        elif variable == 'Precipitation':
            ds  = _load_years('pr')
            da  = (ds['pr'] * 86400.0).clip(min=0.0)   # kg m⁻² s⁻¹ → mm day⁻¹
            da.attrs.update({'units': 'mm/day'})
            out = da.to_dataset(name=variable)

        elif variable == 'ReferenceET':
            # Load the five ET₀ input variables
            ds_tmin = _load_years('tasmin')
            ds_tmax = _load_years('tasmax')
            ds_hurs = _load_years('hurs')
            ds_rsds = _load_years('rsds')
            ds_wind = _load_years('sfcWind')

            tasmin  = ds_tmin['tasmin'] - 273.15        # K → °C
            tasmax  = ds_tmax['tasmax'] - 273.15        # K → °C
            hurs    = ds_hurs['hurs']                   # [%]
            rsds    = ds_rsds['rsds'] * 0.0864          # W m⁻² → MJ m⁻² day⁻¹
            wind10m = ds_wind['sfcWind']                # m s⁻¹ at 10 m

            et0 = _calc_et0_xr(tasmin, tasmax, hurs, rsds, wind10m, elev=elev)
            out = et0.to_dataset(name='ReferenceET')

        else:
            print(f"        WARNING: '{variable}' not recognised – skipping.")
            continue

        # Determine the effective scenario label for the output filename.
        # If the period spans historical and SSP, use the SSP label so users
        # know which projection was used for future years.
        out_scenario = 'historical' if scenario == 'historical' else scenario

        _preproc_and_save(
            out, variable, yearlist, basepath, to_match,
            model, out_scenario, ensemble,
        )
