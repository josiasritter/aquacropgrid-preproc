"""
This script takes user-defined information (domain and time period to be modelled) and downloads and prepares all data needed for running AquaCropGrid
User inputs:
    - path to a vector polygon file (shape) representing the model domain
    - time period to me modelled (start year and end year)
    - desired cell resolution (default 3 arcmin = 0.05 degrees latitude/longitude). Higher resolution may be possible for smaller domains, but keep in mind that the spatial resolution of most input datasets is quite coarse (e.g. 0.25 degrees for future climate data), so higher resolution may not always be useful and may lead to longer processing times and larger file sizes.
    - personal api token from the Copernicus Climate Data Store (https://cds.climate.copernicus.eu/)
Script generates preprocessed datasets (grids) containing:
    - Climate data in daily time step (precipitation, evapotranspiration, minimum temperature, maximum temperature, and initial soil moisture).
      Source is chosen automatically based on the requested time period:
        * Both start_year and end_year are in the past (< current year): AgERA5 reanalysis via Copernicus CDS
        * end_year is in the future (>= current year): NASA NEX-GDDP-CMIP6 climate projections
    - Soil data (content of clay, sand, silt, and soil organic matter) for six soil depth layers from ISCRIC Soilgrids data
    - Crop planting and harvesting calendars (for all supported crop types) from GGCMI data
    - Crop cultivation areas (for all supported crop types) from SPAM data
"""

import os
from preproc_tools import basegrid
from validate_inputs import validate_inputs
# import pdb

## INPUT ARGUMENTS. REPLACE THESE WITH YOUR OWN VALUES
workingdirectory = os.getcwd()   # your home directory
#domain_path = os.path.join(workingdirectory, 'inputdata', 'mekong', 'basin_outline', 'mekong_jrc_outline.geojson')
domain_path = os.path.join(workingdirectory, 'inputdata', 'germany', 'niedersachsen.geojson')   # location and name of your domain shapefile (polygon file representing the model domain). Must be in lat/lon (EPSG:4326) projection.
start_year = 2030
end_year = 2031
cell_resolution = 0.05 # cell resolution in degrees (e.g. 0.05 for 3 arcmin). Resolution of 0.05 degrees is reasonable given the coarse spatial resolution of most input datasets.
api_token = 'xxx'  # your API token when using AgERA5 as climate input, retrieved from your profile page on the Copernicus Climate Data Store (https://cds.climate.copernicus.eu/)

# NASA NEX-GDDP-CMIP6 settings (used for climate projection inputs when end_year >= current year)
nasanex_model    = 'GFDL-CM4'   # CMIP6 model; see https://ds.nccs.nasa.gov/thredds/catalog/AMES/NEX/GDDP-CMIP6/catalog.html
nasanex_scenario = 'ssp245'     # SSP scenario for years >= 2015: 'ssp126', 'ssp245', 'ssp370', 'ssp585'
nasanex_ensemble = 'r1i1p1f1'   # ensemble member (check catalog for model-specific members)

##
def aquacropgrid_preproc(domain_shape_path, start_year, end_year, api_token, cell_resolution=0.05, preprocess=['soil', 'crop_areas', 'cropcalendar', 'climate'],
                        nasanex_model='GFDL-CM4', nasanex_scenario='ssp245', nasanex_ensemble='r1i1p1f1'):
    workingdirectory = os.getcwd()  # your home directory

    # Validate user inputs
    validate_inputs(domain_shape_path, start_year, end_year, api_token)

    # Create template raster file from domain shape for all other datasets to align
    templategrid_path = os.path.join(workingdirectory, 'template_grid.nc')
    to_match, bounds = basegrid(domain_shape_path, cell_resolution, templategrid_path)

    # Download and preprocess soil data from ISRIC Soilgrids
    if 'soil' in preprocess:
        from soil import soil
        soil(domain_shape_path, cell_resolution, workingdirectory, templategrid_path)

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

    # Download and preprocess climate data.
    # Source selection:
    #   AgERA5 reanalysis  – both years within its availability window (1979 to last complete year)
    #   NASA NEX-GDDP-CMIP6 – any other case (start_year < 1979 or end_year >= current year)
    if 'climate' in preprocess:
        import datetime
        current_year = datetime.date.today().year
        AGERA5_START = 1979
        use_agera5 = (start_year >= AGERA5_START) and (end_year < current_year)
        if use_agera5:
            from climate_AgERA5 import climate_AgERA5
            climate_AgERA5(workingdirectory, start_year, end_year, api_token, to_match)
        else:
            from climate_nasanex import climate_nasanex
            climate_nasanex(workingdirectory, start_year, end_year, to_match,
                            model=nasanex_model, scenario=nasanex_scenario, ensemble=nasanex_ensemble)

## Run preprocessing
aquacropgrid_preproc(domain_path, start_year, end_year, api_token, cell_resolution=cell_resolution, preprocess=['soil', 'crop_areas', 'cropcalendar', 'climate'],
                     nasanex_model=nasanex_model, nasanex_scenario=nasanex_scenario, nasanex_ensemble=nasanex_ensemble)
