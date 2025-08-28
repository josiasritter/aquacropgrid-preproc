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

from preproc_tools import preproc_era5, basegrid, makedirs

def climate(basepath, domain_path, start_year, end_year, api_token, cell_resolution, variables=['MinTemp','MaxTemp','Precipitation','ReferenceET','InitSoilwater']):

    ## Years to be downloaded
    yearlist = list(range(start_year, end_year+1))

    # Define area and grid resolution to be downloaded (bounding box)
    templategrid_path = os.path.join(basepath, 'template_grid.nc')
    to_match, bounds = basegrid(domain_path, cell_resolution, templategrid_path) # KEEP HERE FOR NOW IN CASE CLIMATE PREPROCESSING BECOMES SEPARATE PACKAGE
    bounds=[bounds[3],bounds[0],bounds[1],bounds[2]]    # reorder bounds to follow ERA5 CDS definition (N-W-S-E)

    # Prepare download directory
    target_dir = makedirs(basepath, 'rawdata', 'climate')

    # Prepare Copernicus Climate Data Store (CDS) API
    url = 'https://cds.climate.copernicus.eu/api'
    c = cdsapi.Client(url=url, key=api_token)

    ## Download daily min and max temperatures from ERA5-Land daily

    t_stats = ["daily_minimum", "daily_maximum"]
    for stat in t_stats:
        variable = "MinTemp" if stat == "daily_minimum" else "MaxTemp"  # rename variables to AquaCrop definition
        if variable in variables:
            for year in yearlist:   # Split download into separate years
                targetfile = os.path.join(target_dir, variable + str(year) + '.nc')
                if not os.path.exists(targetfile):   # Skip download if file already exists
                    print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(year) + " ***")
                    c.retrieve(
                        'derived-era5-land-daily-statistics',
                        {
                            'variable': ['2m_temperature'],
                            'year': [str(year)],
                            'month': ["01","02","03","04","05","06","07","08","09","10","11","12"],
                            "day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"],
                            "daily_statistic": stat,
                            "time_zone": "utc+00:00",
                            "frequency": "6_hourly",
                            "area": bounds,
                            'data_format': 'netcdf',
                        },
                        targetfile
                    )

            # Combine yearly files into one
            file_paths = [os.path.join(target_dir, variable + str(year) + '.nc') for year in yearlist]
            datasets = [xr.open_dataset(f) for f in file_paths]
            src = xr.concat(datasets, dim='valid_time')
            src = src.sortby('valid_time')    # Make sure all is properly sorted along the time dimension

            # Preprocessing
            preproc_era5(src, variable, yearlist, basepath, to_match)

    ## Download daily precipitation accumulation from ERA5-Land hourly (time step 0000 UTC represents full accumulation of previous day, see https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation)

    variable = 'Precipitation'
    if variable in variables:
        for year in yearlist:
            targetfile = os.path.join(target_dir, variable + str(year) + '.nc')
            if not os.path.exists(targetfile):   # Skip download if file already exists
                print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(year) + " ***")
                c.retrieve(
                    "reanalysis-era5-land",
                    {
                        "variable": ["total_precipitation"],
                        "year": [str(year)],
                        "month": ["01","02","03","04","05","06","07","08","09","10","11","12"],
                        "day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"],
                        "time": ["00:00"],
                        "area": bounds,
                        "data_format": "netcdf",
                        "download_format": "unarchived"
                    },
                    targetfile
                )

        # Last time step (Dec 31) is saved in Jan 1 00:00 of following year -> Separate download necessary with this (annoying) ERA5 definition of accumulation variables in hourly data
        targetfile = os.path.join(target_dir, variable + str(end_year+1) + 'Jan1.nc')
        if not os.path.exists(targetfile):  # Skip download if file already exists
            print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(end_year+1) + "Jan1 ***")
            c.retrieve(
                "reanalysis-era5-land",
                {
                    "variable": ["total_precipitation"],
                    "year": [str(end_year+1)],
                    "month": ["01"],
                    "day": ["01"],
                    "time": ["00:00"],
                    "area": bounds,
                    "data_format": "netcdf",
                    "download_format": "unarchived"
                },
                targetfile
            )

        # Combine yearly files into one
        file_paths = [os.path.join(target_dir, variable + str(year) + '.nc') for year in yearlist] + [os.path.join(target_dir, variable + str(end_year+1) + 'Jan1.nc')]
        datasets = [xr.open_dataset(f) for f in file_paths]
        src = xr.concat(datasets, dim='valid_time')
        src = src.sortby('valid_time')  # Make sure all is properly sorted along the time dimension

        # Preprocessing
        preproc_era5(src, variable, yearlist, basepath, to_match)

    ## Download daily potential evaporation accumulation from ERA5 daily (ERA5 used since ERA5-Land defines potential evaporation differently)
    variable = 'ReferenceET'
    if variable in variables:
        for year in yearlist:
            targetfile = os.path.join(target_dir, variable + str(year) + '.nc')
            if not os.path.exists(targetfile):  # Skip download if file already exists
                print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(year) + " ***")
                c.retrieve(
                    "derived-era5-single-levels-daily-statistics",
                    {
                        "product_type": "reanalysis",
                        "variable": ["potential_evaporation"],
                        "year": [str(year)],
                        "month": ["01","02","03","04","05","06","07","08","09","10","11","12"],
                        "day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"],
                        "daily_statistic": "daily_sum",
                        "time_zone": "utc+00:00",
                        "frequency": "1_hourly",
                        "area": bounds,
                        "data_format": "netcdf"
                    },
                    targetfile
                )

        # Combine yearly files into one
        file_paths = [os.path.join(target_dir, variable + str(year) + '.nc') for year in yearlist]
        datasets = [xr.open_dataset(f) for f in file_paths]
        src = xr.concat(datasets, dim='valid_time')
        src = src.sortby('valid_time')    # Make sure all is properly sorted along the time dimension

        # Preprocessing
        preproc_era5(src, variable, yearlist, basepath, to_match)

    ## Download Volumetric soil water content [m3/m3] for initial time step (has four soil depth layers) from ERA5-Land hourly

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
        preproc_era5(src, variable, yearlist, basepath, to_match)
