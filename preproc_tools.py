# IMPORT LIBRARIES

import xarray as xr
from rasterio.warp import Resampling
import numpy as np
#from scipy.ndimage import convolve
import pandas as pd
import geopandas as gpd
import os
import glob
import time
import requests
from zipfile import ZipFile
import pdb # pdb.set_trace()

import rioxarray
from shapely.geometry import mapping, box

## Download functions

def download_url(url, download_path="./"):
    """Download data from a given url to a given download path on the hard drive"""
    while True:
        try:
            response = requests.get(url, timeout=30)
        except:
            print("    *** WARNING! Download froze, likely due to network instability. Restarting download...")
            time.sleep(2)
        else:
            if response.ok:
                # inmemory = tiff.imread(io.BytesIO(response.content))  # reads download data directly to memory (without writing to disk)
                with open(download_path, 'wb') as savefile:
                    savefile.write(response.content)  # write download data to disk
                break

def unzip_all(dir='.'):
    """Iteratively unzip files in given directory and subdirectories until no zip files remain"""
    zipfiles = True
    while bool(zipfiles):
        zipfiles = []
        for root, d_names, f_names in os.walk(dir):
            for f in f_names:
                if f.endswith('.zip'):
                    zipfiles.append(os.path.join(root, f))
        if bool(zipfiles):
            for file in zipfiles:
                targetdir = file[:-4]
                if not os.path.exists(targetdir):
                    os.mkdir(targetdir)
                with ZipFile(file, 'r') as zObject:
                    zObject.extractall(path=targetdir)
                os.remove(file)

## Preprocessing functions
def preproc_spam(basepath, download_dir, refyear, spam_variable, domain_path, to_match):
    mask = gpd.read_file(domain_path)
    to_match.rio.write_crs(4326, inplace=True)

    # Dictionary that connects FAO crop ID's to crop names in AquaCrop. Used for layer naming
    crop_dict = {'BARL': 'Barley','COTT': 'Cotton','BEAN': 'DryBean','MAIZ': 'Maize','RICE': 'PaddyRice','POTA': 'Potato','SORG': 'Sorghum','SOYB': 'Soybean','SUGB': 'SugarBeet','SUGC': 'SugarCane','SUNF': 'Sunflower','WHEA': 'Wheat_summer','CASS': 'Cassava'}

    # Dictionnary to rename technology
    tech_dict = {'R': 'rf', 'I': 'ir'}

    # CHECK IF ALL FOUR VARIABLES ARE NEEDED. PROBABLY JUST PHYSICAL AREA AND - IF CALBRATION INCLUDED - YIELD
    search_criteria = os.path.join(download_dir, '**', '*.tif')
    layers = glob.glob(search_criteria, recursive=True)
    file_to_mosaic = []

    for layer in layers:

        # Preserve only layers referring to rainfed (R) or irrigated (I) technology
        technique = layer[-5:-4]
        if technique not in ['R','I']:
            os.remove(layer)
            continue

        # Get crop ID from layer name and preserve only if crop type is supported by AquaCrop
        cropID = layer[-10:-6]
        if cropID not in crop_dict:
            os.remove(layer)
            continue

        print(cropID, technique)

        src = xr.load_dataset(layer)
        src_proj = src.rio.reproject_match(to_match, resampling=Resampling.average) # changed from nearest-neighbor to area-weighted resampling to reduce distortions 
        #src_clip = src_proj.rio.clip(mask.geometry.apply(mapping))
        src_clip = safe_clip(src_proj, mask)

        # Rename variables from FAO crop ID to crop strings used by AquaCrop
        layername_base = crop_dict.get(cropID) + '_' + tech_dict.get(technique)

        # Change raster values due to resolution change from 5 arcmin to 3 arcmin for physical and harvested area and production rasters
        if (spam_variable == 'physical_area') | (spam_variable == 'harvested_area') | (spam_variable == 'production'):
            src_clip[layername_base + '_' + spam_variable] = src_clip['band_data'] * 3**2 / 5**2 # TESTED FOR MEKONG: OVERALL CROP AREA FOR IRRIGATED RICE DECREASES BY ABOUT 1% -> negligible!
            #src_clip[layer[-12: -4] + '_percentage'] = src_clip[layer[-12: -4] + '_area'] / meters.ha  # Percentage only needed to restrict model to cells with very small crop area

        # Renaming variable for yield and change from kg/ha to t/ha
        elif spam_variable == 'yield':
            src_clip[layername_base + '_yield'] = src_clip['band_data'] / 1000

        # Drop band_data and append file
        src_clip = src_clip.drop_vars(['band_data'])
        file_to_mosaic.append(src_clip)

    # Merge data into one file
    src_mosaic = xr.merge(file_to_mosaic)
    target_dir = makedirs(basepath, 'processed', '')
    targetfile = os.path.join(target_dir, 'spam' + refyear + '_' + spam_variable + '.nc')
    src_mosaic.to_netcdf(targetfile)

def spam_refyear(start_year, end_year):
    # Choose SPAM data reference year (available reference years are 2010 or 2020) according to average of start and end year of modelling horizon
    import numpy as np
    avg_year = np.mean([start_year, end_year])
    avg_year = np.ceil(avg_year)
    SPAM_refyears = [2010, 2020]
    refyear = min(SPAM_refyears, key=lambda x: abs(x - avg_year)) # select nearest reference year available in SPAM
    refyear = str(refyear)

    return refyear

## Preprocessing climate data
def preproc_era5(src, variable, yearlist, basepath, to_match):
    print("        *** PREPROCESSING CLIMATE DATA: " + variable + " ***")

    # Preparations
    src = ensure_xy_dims(src)
    src.rio.write_crs(4326, inplace=True)

    varname_dict = {'MinTemp': 't2m', 'MaxTemp': 't2m', 'Precipitation': 'tp', 'ReferenceET': 'pev'}  # Names of data variables in AquaCrop and ERA5 data, respectively

    # Adjust units and data variable names to AquaCrop definitions
    if variable in ['MinTemp', 'MaxTemp']:
        src = src.rename_vars({varname_dict.get(variable): variable})  # Rename data variables to AquaCrop names
        src = src - 273.15  # Convert from Kelvin to Celsius
    if variable in ['Precipitation', 'ReferenceET']:
        src = src.rename_vars({varname_dict.get(variable): variable})  # Rename data variables to AquaCrop names
        src = src * 1000 # Convert from metres to millimetres
        if variable in ['ReferenceET']:
            src = -src  # Multiply by -1 since raw data is negative due to ERA5 definition (i.e. downward fluxes are positive, upward fluxes negative)
        elif variable in ['Precipitation']:
            # Adjust time coordinate to date of previous day, as Jan 1 00:00 represents daily accumulation of Dec 31 of previous year (see ERA5 time step definition: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation )
            src = src.isel(valid_time=slice(1, None))   # Delete first time step as it represents daily accumulation of 31/12 of previous year
            valid_dates = pd.to_datetime(src.valid_time.values)  # Convert valid_time to pandas DatetimeIndex
            new_dates = (valid_dates.normalize() - pd.Timedelta(days=1))  # Strip time part (keep just the date) and shift back by one day
            src = src.assign_coords(valid_time=new_dates)  # Assign adjusted timestamps back to the dataset
            src = src.drop_vars('expver')
        src[variable] = src[variable].where(src[variable] >= 0, 0)  # Set negative precipitation and evaporation values to 0.
    if variable in ['InitSoilwater']:
        src = src.drop_vars(['expver','number'])

    # Fill missing values in data variable(s)
    import scipy.ndimage as ndi
    variables = list(src.data_vars)
    for variab in variables:    # Loop to make it work also for InitSoilwater data, which has four data variables (corresponding to four soil depths)
        data3d = src[variab].to_numpy()
        nodata_mask = np.isnan(data3d[0])  # This assumes that all timesteps have the same missing cells (fair assumption, since the same ERA5-land water mask is applied to all time steps)
        dist, nearest_indices = ndi.distance_transform_edt(nodata_mask, return_indices=True)  # Nearest-neighbour interpolation
        src[variab].data = data3d[:, nearest_indices[0], nearest_indices[1]]

    # Resample to project grid and mask
    src_reproj = src.rio.reproject_match(to_match, resampling=Resampling.nearest)
    src_reproj = src_reproj.rename({'valid_time': 'time', 'x': 'longitude', 'y': 'latitude'})
    src_masked = src_reproj.where(to_match['Band1'] == 1)

    #if variable in ['InitSoilwater']:  # NOT NEEDED AFTER ALL BECAUSE AQUACROP CAN DO IT USING ORIGINAL SOIL LAYERS FROM ERA5 (see function get_initial_WC in aquacropgrid.py)
    #    src_reproj = convert_soildepthlayers(src_reproj)

    # Prepare output directory
    target_dir = makedirs(basepath, 'processed', '')
    targetfile = os.path.join(target_dir, variable + str(yearlist[0]) + str(yearlist[-1]) + '.nc')

    # Save to disk
    src_masked.to_netcdf(targetfile, mode='w', encoding = {variables[0]: {'zlib': True, 'complevel': 4}})  # Save to disk with moderate compression


## Helper functions

def basegrid(domain_shape_path, resolution, templategrid_path):    # Creates basic raster file as template for all preprocessing scripts (i.e. what is read in all other scripts as "to_match" file)
    from rasterio.transform import from_origin
    from rasterio import features
    from affine import Affine

    # Read shapefile
    mask = gpd.read_file(domain_shape_path)

    # Check if all geometries are Polygon or MultiPolygon
    if not all(g in ['Polygon', 'MultiPolygon'] for g in mask.geom_type.unique()):
        raise Exception("The input polygon file contains geometries other than Polygon/MultiPolygon.")

    # Check if it's EPSG:4326
    if not mask.crs.to_epsg() == 4326:
        raise Exception("The polygon file must be in EPSG:4326.")

    full_geom = mask.unary_union
    xmin, ymin, xmax, ymax = full_geom.bounds # in lat/lon
    xmin, ymin = np.floor(xmin * (1/resolution)) / (1/resolution) , np.floor(ymin * (1/resolution)) / (1/resolution) # round bounds down to cell resolution
    xmax, ymax = np.ceil(xmax * (1/resolution)) / (1/resolution) , np.ceil(ymax * (1/resolution)) / (1/resolution) # round upper bounds up to cell resolution
    bounds = [xmin, ymin, xmax, ymax]

    if os.path.exists(templategrid_path):
        ds = xr.open_dataset(templategrid_path)

    else:
        w = round((xmax - xmin) / resolution)
        h = round((ymax - ymin) / resolution)
        lons = np.arange(xmin+resolution/2, xmax, resolution)
        lats = np.arange(ymin+resolution/2, ymax, resolution)

        # Affine transform: top-left origin
        transform = Affine.translation(xmin, ymax) * Affine.scale(resolution, -resolution)

        # Mask: True = inside polygon
        mask = features.geometry_mask([full_geom],out_shape=(h, w),transform=transform,invert=True)

        # Create data array: inside = 1, outside = nodata
        nodata_value = 0
        data = np.full((h, w), nodata_value, dtype=np.uint8)
        data[mask] = 1

        # CF convention prefers descending latitudes
        lats = lats[::-1]

        # Create xarray Dataset
        ds = xr.Dataset({"Band1": (["latitude", "longitude"], data),},coords={"latitude": lats,"longitude": lons,},attrs={"Conventions": "CF-1.8","title": "Template raster from polygon shapefile","crs": "EPSG:4326"})
        ds["spatial_ref"] = xr.DataArray(0, attrs={"grid_mapping_name": "latitude_longitude", "epsg_code": 4326, "semi_major_axis": 6378137.0, "inverse_flattening": 298.257223563, "long_name": "CRS definition"})
        ds["Band1"].attrs.update({"grid_mapping": "spatial_ref", "_FillValue": nodata_value, "missing_value": nodata_value})

        # Save to NetCDF
        ds.to_netcdf(templategrid_path)

    return ds, bounds

# Ensures the dataset has 'x' and 'y' dimensions (for resampling in xarray)
def ensure_xy_dims(ds):
    # Rename dimensions to 'x' and 'y' if they are named differently
    if 'longitude' in ds.dims and 'latitude' in ds.dims:
        ds = ds.rename({'longitude': 'x', 'latitude': 'y'})
    elif 'lon' in ds.dims and 'lat' in ds.dims:
        ds = ds.rename({'lon': 'x', 'lat': 'y'})
    return ds

def makedirs(basepath, level1, level2):
    level1_dir = os.path.join(basepath, level1)
    if not os.path.exists(level1_dir):
        os.mkdir(level1_dir)
    level2_dir = os.path.join(level1_dir, level2)
    if not os.path.exists(level2_dir):
        os.mkdir(level2_dir)
    return level2_dir


def safe_clip(src, mask):
    """
    Clip a raster (src) using geometries from mask, skipping polygons
    that fall entirely in nodata regions or outside the raster extent.

    Parameters
    ----------
    src : rioxarray.DataArray
        The source raster opened with rioxarray.
    mask : geopandas.GeoDataFrame
        Polygon(s) to clip with. Can be MultiPolygon.

    Returns
    -------
    xarray.DataArray
        Clipped raster with same CRS and structure as src.
    """
    # Ensure CRS match
    if mask.crs != src.rio.crs:
        mask = mask.to_crs(src.rio.crs)

    # Explode multipolygons
    mask_exploded = mask.explode(ignore_index=True)

    # Get raster bounds as shapely box
    raster_bounds = box(*src.rio.bounds())

    valid_clips = []

    for i, geom in enumerate(mask_exploded.geometry):
        # Skip polygons completely outside raster extent
        if not geom.intersects(raster_bounds):
            continue
        try:
            clipped = src.rio.clip([mapping(geom)], src.rio.crs, drop=True)
            # Only keep if it actually contains data
            if clipped.notnull().any():
                valid_clips.append(clipped)
        except rioxarray.exceptions.NoDataInBounds:
            continue

    # If no valid clips, return empty raster with same structure
    if not valid_clips:
        print("Warning: no valid polygons contained data.")
        empty = src.copy(deep=True)
        empty[:] = src.rio.nodata
        return empty

    # Merge the valid clipped parts into one raster
    combined = xr.concat(valid_clips, dim="band").max(dim="band")

    # Preserve metadata & CRS
    combined.rio.write_crs(src.rio.crs, inplace=True)
    combined.attrs.update(src.attrs)

    return combined