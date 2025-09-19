"""
This script takes user-defined information (domain and time period to be modelled) and downloads and prepares all data needed for running AquaCropGrid
User inputs:
    - path to a vector polygon file (shape) representing the model domain
    - time period to me modelled (start year and end year)
    - personal api token from the Copernicus Climate Data Store (https://cds.climate.copernicus.eu/)
    - desired cell resolution (default 3 arcmin = 0.05 degrees latitude/longitude). For now, only 0.05 degrees works. User-defined cell size will be added in the future.
Script generates preprocessed datasets (grids) containing:
    - Climate data in daily time step (precipitation, evapotranspiration, minimum temperature, maximum temperature, and initial soil moisture) from ERA5 and ERA5-land data
    - Soil data (content of clay, sand, silt, and soil organic matter) for six soil depth layers from ISCRIC Soilgrids data
    - Crop cultivation areas (for all supported crop types) from SPAM data
    - Crop planting and harvesting calendars (for all supported crop types) from GGCMI data
"""

import os
from preproc_tools import basegrid

## INPUT ARGUMENTS. REPLACE THESE WITH YOUR OWN VALUES
workingdirectory = os.getcwd()   # your home directory
domain_path = os.path.join(workingdirectory, 'inputdata', 'mekong', 'basin_outline', 'mekong_jrc_outline.geojson')
start_year = 2014
end_year = 2015
api_token = 'xx'  # your API token, retrieved from your profile page on the Copernicus Climate Data Store (https://cds.climate.copernicus.eu/)

##
def aquacropgrid_preproc(domain_shape_path, start_year, end_year, api_token, cell_resolution=0.05, preprocess=['soil', 'crop_areas', 'cropcalendar', 'climate']):
    workingdirectory = os.getcwd()  # your home directory

    # Creat template raster file from domain shape for all other datasets to align
    templategrid_path = os.path.join(workingdirectory, 'template_grid.nc')
    to_match, bounds = basegrid(domain_path, cell_resolution, templategrid_path)


    # Download and preprocess soil data from ISRIC Soilgrids
    if 'soil' in preprocess:
        from soil import soil
        soil(domain_shape_path, cell_resolution, workingdirectory)

    # Download and preprocess crop areas (crop mask) and crop yield from SPAM data (https://www.mapspam.info/)
    if 'crop_areas' in preprocess:
        from crop_areas import crop_areas
        spam_variable = 'physical_area' # crop masks, seperately for rainfed and irrigated areas
        crop_areas(domain_shape_path, spam_variable, start_year, end_year, workingdirectory, to_match)
        # spam_variable = 'yield' # crop yields, seperately for rainfed and irrigated areas. Used only for calibration and/or validation
        # crop_areas(domain_shape_path, spam_variable, start_year, end_year, workingdirectory, to_match)

    # Download and preprocess crop calendar from GGCMI (https://zenodo.org/records/5062513)
    if 'cropcalendar' in preprocess:
        from cropcalendar_module import cropcalendar 
        cropcalendar(domain_shape_path, workingdirectory, templategrid_path)

    # Download and preprocess climate data and initial soil moisture from ERA5 and ERA5-Land
    if 'climate' in preprocess:
        from climate import climate
        climate(workingdirectory, domain_shape_path, start_year, end_year, api_token, cell_resolution)

## Run preprocessing
aquacropgrid_preproc(domain_path, start_year, end_year, api_token, preprocess=['soil', 'crop_areas', 'cropcalendar', 'climate'])
