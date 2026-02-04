"""
This script downloads and prepocesses the following data from ERA5 and ERA5-Land:
    - Four daily climate variables:
        - Minimum temperature [K]
        - Maximum temperature [K]
        - Precipitation [m]
        - Potential evapotranspiration [m]
    - Soil moisture of initial time step [m3/m3]
"""

import xarray as xr
#import rioxarray as rio
#import numpy as np
#from rasterio.warp import Resampling
import os
import geopandas as gpd
import cdsapi
from shapely.geometry import mapping
import pdb

from preproc_tools import agera5_merge_yearly, preproc_era5, preproc_agera5, basegrid, makedirs, unzip_all

def climate_AgERA5(basepath, domain_path, start_year, end_year, api_token, cell_resolution, variables=['MinTemp','MaxTemp','Precipitation','ReferenceET','InitSoilwater']):

    ## Years to be downloaded
    yearlist = list(range(start_year, end_year+1))

    # Define area and grid resolution to be downloaded (bounding box)
    templategrid_path = os.path.join(basepath, 'template_grid.nc')
    to_match, bounds = basegrid(domain_path, cell_resolution, templategrid_path) # KEEP HERE FOR NOW IN CASE CLIMATE PREPROCESSING BECOMES SEPARATE PACKAGE
    bounds=[bounds[3],bounds[0],bounds[1],bounds[2]]    # reorder bounds to follow ERA5 CDS definition (N-W-S-E)

    # Prepare download directory
    target_dir = makedirs(basepath, 'rawdata', 'climate')

    # Prepare variable names and stats for API request
    varname_api = {'MinTemp': '2m_temperature', 'MaxTemp': '2m_temperature', 'Precipitation': 'precipitation_flux', 'ReferenceET': 'reference_evapotranspiration'}  # Names of data variables in AquaCrop and AgERA5 api, respectively
    stats_api = {'MinTemp': '24_hour_minimum', 'MaxTemp': '24_hour_maximum', 'Precipitation': '', 'ReferenceET': ''}  # Stats to be requested from API for each variable. Only needed for min and max temperature, as AgERA5 api provides daily accumulations for precipitation and reference ET

    # Prepare Copernicus Climate Data Store (CDS) API
    url = 'https://cds.climate.copernicus.eu/api'
    c = cdsapi.Client(url=url, key=api_token)

    ## Download Mintemp, Maxtemp, ReferenceET, and Precipitation in daily timestep from AgERA5.
    for variable in variables:
        if variable == 'InitSoilwater': # Skip initial soil water here, as it is downloaded from ERA5-Land below
            continue
        for year in yearlist:   # Split downloads into yearly chunks to avoid large files
            targetfile = os.path.join(target_dir, variable + str(year) + '.zip')
            yearfile = os.path.join(target_dir, variable + str(year) + '.nc')
            if not os.path.exists(targetfile) and not os.path.exists(yearfile):  # Skip download if zip file or merged yearly .nc file already exist
                print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(year) + " ***")
                c.retrieve(
                    "sis-agrometeorological-indicators",
                    {
                        "variable": [varname_api.get(variable)],
                        "statistic": [stats_api.get(variable)],
                        "year": [str(year)],
                        "month": ["01","02","03","04","05","06","07","08","09","10","11","12"],
                        "day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"],
                        "area": bounds,
                        "version": "2_0",
                    },
                    targetfile
                )

            # Unzip download file and merge the resulting daily files into one yearly file (AgERA5 api returns a zip file containing daily .nc files)
            if not os.path.exists(yearfile):  # Skip unzipping and file merging if yearly .nc file already exists
                agera5_merge_yearly(target_dir, yearfile)

        # Combine yearly files into one
        file_paths = [os.path.join(target_dir, variable + str(year) + '.nc') for year in yearlist]
        datasets = [xr.open_dataset(f) for f in file_paths]
        src = xr.concat(datasets, dim='time')
        src = src.sortby('time')    # Make sure all is properly sorted along the time dimension

        # Preprocessing
        preproc_agera5(src, variable, yearlist, basepath, to_match)


    ## Download Volumetric soil water content [m3/m3] for initial time step (has four soil depth layers) from ERA5-Land hourly (not available in AgERA5)
    variable = 'InitSoilwater'
    if variable in variables:
        targetfile = os.path.join(target_dir, variable + str(start_year) + '.nc')
        if not os.path.exists(targetfile):  # Skip download if file already exists
            print("        *** DOWNLOADING CLIMATE DATA: " + variable + " ***")
            c.retrieve(
                "reanalysis-era5-land",
                {
                    "variable": [
                        "volumetric_soil_water_layer_1",
                        "volumetric_soil_water_layer_2",
                        "volumetric_soil_water_layer_3",
                        "volumetric_soil_water_layer_4"
                    ],
                    "year": [str(start_year)],
                    "month": ["01"],
                    "day": ["01"],
                    "time": ["00:00"],
                    "data_format": "netcdf",
                    "download_format": "unarchived",
                    "area": bounds
                },
                targetfile
            )

        # Preprocessing
        src = xr.open_dataset(targetfile)
        #preproc_era5(src, variable, yearlist, basepath, to_match)
        preproc_agera5(src, variable, yearlist, basepath, to_match)
