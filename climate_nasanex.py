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
    # Cast all inputs to float32 so intermediates stay in float32, halving memory vs float64
    f = np.float32
    tasmin  = tasmin.astype(f)
    tasmax  = tasmax.astype(f)
    hurs    = hurs.astype(f)
    rsds    = rsds.astype(f)
    wind10m = wind10m.astype(f)

    # --- Physical constants (float32 to keep arithmetic in float32) ----------
    Gsc  = f(0.0820)    # solar constant              [MJ m⁻² min⁻¹]
    sbc  = f(4.903e-9)  # Stefan-Boltzmann constant   [MJ m⁻² day⁻¹ K⁻⁴]
    alb  = f(0.23)      # grass albedo (FAO-56)
    Cn   = f(900.0)     # short-crop numerator constant
    Cd   = f(0.34)      # short-crop denominator constant
    G    = f(0.0)       # soil heat flux ≈ 0 for daily time step

    # --- Atmospheric pressure and psychrometric constant ---------------------
    AtmP = f(101.3) * ((f(293.0) - f(0.0065) * f(elev)) / f(293.0)) ** f(5.26)  # [kPa]
    psy  = f(0.000665) * AtmP                                                      # [kPa °C⁻¹]

    # --- Wind speed: 10 m → 2 m  (FAO-56 Eq. 47) ----------------------------
    u2 = wind10m * f(4.87 / np.log(67.8 * 10.0 - 5.42))  # pre-compute scalar as float32

    # --- Temperature and humidity --------------------------------------------
    tmean = (tasmax + tasmin) * f(0.5)
    rh    = hurs.clip(f(0.0), f(100.0))

    # Saturation vapour pressure [kPa]
    e0max = f(0.6108) * np.exp(f(17.27) * tasmax / (tasmax + f(237.3)))
    e0min = f(0.6108) * np.exp(f(17.27) * tasmin / (tasmin + f(237.3)))
    es    = (e0max + e0min) * f(0.5)

    # Actual vapour pressure [kPa]
    ea = (rh * f(0.01)) * es

    # Slope of saturation vapour pressure curve [kPa °C⁻¹]
    Delta = (
        f(4098.0) * f(0.6108) * np.exp(f(17.27) * tmean / (tmean + f(237.3)))
    ) / (tmean + f(237.3)) ** f(2.0)

    # --- Radiation -----------------------------------------------------------
    # Julian day of year (1-D along 'time')
    J = tasmin.time.dt.dayofyear.astype(np.float32)

    # Latitude in radians (1-D along 'y')
    lat_rad = np.deg2rad(tasmin.y.astype(np.float32))

    # Inverse relative Earth-Sun distance and solar declination
    pi365    = f(2.0 * np.pi / 365.0)  # pre-computed float32 constant
    dr       = f(1.0) + f(0.033) * np.cos(pi365 * J)
    sol_decl = f(0.409) * np.sin(pi365 * J - f(1.39))

    # Sunset hour angle [rad]:  (-tan φ · tan δ) must be clamped to [-1, 1]
    # xarray broadcasts (y,) × (time,)  →  (time, y) automatically
    arccos_arg = (-np.tan(lat_rad) * np.tan(sol_decl)).clip(f(-1.0), f(1.0))
    ws = np.arccos(arccos_arg)  # shape (time, y)

    # Extraterrestrial radiation [MJ m⁻² day⁻¹]  shape (time, y)
    Ra = (
        f(24.0 * 60.0 / np.pi) * Gsc * dr * (
            ws * np.sin(lat_rad) * np.sin(sol_decl)
            + np.cos(lat_rad) * np.cos(sol_decl) * np.sin(ws)
        )
    )

    # Net shortwave radiation [MJ m⁻² day⁻¹]  shape (time, y, x)
    Rns = (f(1.0) - alb) * rsds

    # Clear-sky radiation [MJ m⁻² day⁻¹];  avoid division by zero (polar night)
    Rs0 = (f(0.75) + f(2.0e-5) * f(elev)) * Ra
    Rs0 = Rs0.where(Rs0 > f(0.0), other=f(np.nan))

    # Cloudiness fraction, clamped to [0.05, 1.0]
    fcd = (f(1.35) * rsds / Rs0 - f(0.35)).clip(f(0.05), f(1.0))
    fcd = fcd.where(Ra > f(0.0), other=f(0.05))  # set to minimum during polar night

    # Net longwave radiation [MJ m⁻² day⁻¹]
    Rnl = (
        sbc * fcd
        * (f(0.34) - f(0.14) * np.sqrt(ea.clip(f(0.0))))
        * (((tasmax + f(273.16)) ** f(4.0) + (tasmin + f(273.16)) ** f(4.0)) * f(0.5))
    )

    # Net radiation [MJ m⁻² day⁻¹]
    Rn = Rns - Rnl

    # --- FAO-56 PM ET₀ [mm day⁻¹] -------------------------------------------
    et0 = (
        f(0.408) * Delta * (Rn - G)
        + psy * (Cn / (tmean + f(273.0))) * u2 * (es - ea)
    ) / (Delta + psy * (f(1.0) + Cd * u2))

    # Set ET₀ to 0 during polar night and floor at 0
    et0 = et0.where(Ra > f(0.0), other=f(0.0)).clip(f(0.0))

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
        files = [
            os.path.join(target_dir, f"{nex_var}_{model}_{_scenario_for_year(year, scenario)}_{year}.nc")
            for year in yearlist
        ]
        ds = xr.open_mfdataset(files, combine='by_coords').sortby('time')
        # Normalise any 0-360 longitudes to -180/180
        ds = ensure_xy_dims(ds)
        if 'x' in ds.coords and float(ds.x.max()) > 180.0:
            ds = ds.assign_coords(x=((ds.x + 180.0) % 360.0 - 180.0)).sortby('x')
        return ds

    def _load_one_year(nex_var, year):
        """Open a single year file lazily (used for year-by-year ET₀ to limit peak RAM)."""
        scen = _scenario_for_year(year, scenario)
        path = os.path.join(target_dir, f"{nex_var}_{model}_{scen}_{year}.nc")
        ds = xr.open_dataset(path)
        ds = ensure_xy_dims(ds)
        if 'x' in ds.coords and float(ds.x.max()) > 180.0:
            ds = ds.assign_coords(x=((ds.x + 180.0) % 360.0 - 180.0)).sortby('x')
        return ds

    # ------------------------------------------------------------------
    # Step 3 – Process each requested output variable.
    # All data is cast to float32 (raw NEX files are float32; Python float
    # arithmetic would otherwise silently upcast to float64).
    # ReferenceET is computed year-by-year to cap peak RAM: loading all 5
    # input variables × all years at once is prohibitive for large domains.
    # ------------------------------------------------------------------
    out_scenario = 'historical' if scenario == 'historical' else scenario

    for variable in variables:

        if variable in ('MinTemp', 'MaxTemp'):
            nex_var = _direct[variable]
            ds = _load_years(nex_var)
            # np.float32 constant keeps subtraction in float32
            da = ds[nex_var].astype(np.float32) - np.float32(273.15)
            da.attrs.update({'units': 'degC'})
            out = da.to_dataset(name=variable)
            _preproc_and_save(out, variable, yearlist, basepath, to_match, model, out_scenario, ensemble)
            ds.close()

        elif variable == 'Precipitation':
            ds = _load_years('pr')
            da = (ds['pr'].astype(np.float32) * np.float32(86400.0)).clip(min=np.float32(0.0))
            da.attrs.update({'units': 'mm/day'})
            out = da.to_dataset(name=variable)
            _preproc_and_save(out, variable, yearlist, basepath, to_match, model, out_scenario, ensemble)
            ds.close()

        elif variable == 'ReferenceET':
            # Year-by-year: only 5 × 1-year float32 arrays (+intermediates) in RAM at once
            print("        Computing ET₀ year-by-year to limit peak memory...")
            et0_parts = []
            for year in yearlist:
                ds_tmin = _load_one_year('tasmin', year)
                ds_tmax = _load_one_year('tasmax', year)
                ds_hurs = _load_one_year('hurs',   year)
                ds_rsds = _load_one_year('rsds',   year)
                ds_wind = _load_one_year('sfcWind', year)

                et0_yr = _calc_et0_xr(
                    ds_tmin['tasmin'].astype(np.float32) - np.float32(273.15),
                    ds_tmax['tasmax'].astype(np.float32) - np.float32(273.15),
                    ds_hurs['hurs'].astype(np.float32),
                    ds_rsds['rsds'].astype(np.float32) * np.float32(0.0864),
                    ds_wind['sfcWind'].astype(np.float32),
                    elev=elev,
                )
                # Materialise and immediately release input file handles
                et0_parts.append(et0_yr.load().astype(np.float32))
                for _ds in [ds_tmin, ds_tmax, ds_hurs, ds_rsds, ds_wind]:
                    _ds.close()

            et0 = xr.concat(et0_parts, dim='time').sortby('time')
            out = et0.to_dataset(name='ReferenceET')
            _preproc_and_save(out, variable, yearlist, basepath, to_match, model, out_scenario, ensemble)

        else:
            print(f"        WARNING: '{variable}' not recognised – skipping.")
