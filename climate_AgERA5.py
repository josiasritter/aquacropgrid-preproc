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
#import pdb # pdb.set_trace()

from preproc_tools import agera5_merge_yearly, preproc_agera5, basegrid, makedirs, unzip_all

import socket
import requests
import dns.resolver  # pip install dnspython

# Set of functions to automatically overcome DNS-based errors, often found on university networks
def force_resolve(ip, hostname="cds.climate.copernicus.eu"):
    """
    Force a specific hostname to resolve to a given IP address
    inside Python, without touching system DNS.
    """
    print(f'Forcing resolve: {hostname} -> {ip}')
    orig_getaddrinfo = socket.getaddrinfo
 
    def new_getaddrinfo(*args, **kwargs):
        if args[0] == hostname:
            return orig_getaddrinfo(ip, *args[1:], **kwargs)
        return orig_getaddrinfo(*args, **kwargs)
 
    socket.getaddrinfo = new_getaddrinfo
 
 
def resolve_via_doh(hostname="cds.climate.copernicus.eu", doh_url="https://cloudflare-dns.com/dns-query"):
    """
    Resolve hostname via DNS-over-HTTPS, bypassing any
    port-53 interception by the university network.
    """
    resp = requests.get(
        doh_url,
        params={"name": hostname, "type": "A"},
        headers={"Accept": "application/dns-json"},
    )
    resp.raise_for_status()
    data = resp.json()
    # Filter for A records (type 1)
    a_records = [ans["data"] for ans in data.get("Answer", []) if ans["type"] == 1]
    if not a_records:
        raise RuntimeError(f"No A records found for {hostname} via DoH")
    ip = a_records[0]
    print(f"Resolved {hostname} to {ip} via DoH")
    return ip
 
 
_dns_fallback_applied = False  # Module-level flag so we only apply the monkey-patch once
 
 
def _is_dns_error(exc):
    """Check whether an exception (or its chain) is a DNS resolution failure."""
    # Walk the cause chain — cdsapi wraps errors in requests exceptions
    current = exc
    while current is not None:
        if isinstance(current, socket.gaierror):
            return True
        # requests wraps socket errors in ConnectionError
        if 'Name or service not known' in str(current) or 'getaddrinfo failed' in str(current):
            return True
        current = getattr(current, '__cause__', None) or getattr(current, '__context__', None)
    return False
 
 
def retrieve_with_dns_fallback(client, dataset, request, target):
    """
    Wrapper around cdsapi retrieve that automatically falls back to
    public DNS resolution if the university network can't resolve the CDS hostname.
    """
    global _dns_fallback_applied
 
    # If we've already applied the fix in a previous call, just go straight through
    if _dns_fallback_applied:
        return client.retrieve(dataset, request, target)
 
    try:
        return client.retrieve(dataset, request, target)
    except Exception as e:
        if _is_dns_error(e):
            print('\nDNS resolution failed — falling back to public DNS...')
            ip = resolve_via_doh()
            force_resolve(ip)
            _dns_fallback_applied = True
            return client.retrieve(dataset, request, target)
        else:
            raise  # Not a DNS problem — re-raise as-is
            
def ensure_cds_dns(hostname="cds.climate.copernicus.eu"):
    """
    Check if the CDS hostname resolves via system DNS.
    If not, fall back to public DNS and monkey-patch.
    """
    try:
        socket.getaddrinfo(hostname, 443)
        print(f'DNS resolution OK for {hostname}')
    except socket.gaierror:
        print(f'System DNS failed for {hostname} — falling back to public DNS...')
        ip = resolve_via_doh(hostname)
        force_resolve(ip, hostname)


# Continue with main script functionality
def climate_AgERA5(basepath, start_year, end_year, api_token, to_match, variables=['MinTemp','MaxTemp','Precipitation','ReferenceET','InitSoilwater']):

    ## Years to be downloaded
    yearlist = list(range(start_year, end_year+1))

    # Define area and grid resolution to be downloaded (bounding box)
    templategrid_path = os.path.join(basepath, 'template_grid.nc')
    _tpl = xr.open_dataset(templategrid_path)           # Read spatial extent from template grid file
    _tpl.rio.write_crs(4326, inplace=True)
    bounds = list(_tpl.rio.bounds())                    # [xmin, ymin, xmax, ymax]
    bounds = [round(b,2) for b in bounds]               # round coordinates to shorten filenames (Windows limitation)
    bounds=[bounds[3],bounds[0],bounds[1],bounds[2]]    # reorder bounds to follow ERA5 CDS definition (N-W-S-E)

    # Prepare download directory
    target_dir = makedirs(basepath, 'rawdata', 'climate')

    # Prepare variable names and stats for API request
    varname_api = {'MinTemp': '2m_temperature', 'MaxTemp': '2m_temperature', 'Precipitation': 'precipitation_flux', 'ReferenceET': 'reference_evapotranspiration'}  # Names of data variables in AquaCrop and AgERA5 api, respectively
    stats_api = {'MinTemp': '24_hour_minimum', 'MaxTemp': '24_hour_maximum', 'Precipitation': '', 'ReferenceET': ''}  # Stats to be requested from API for each variable. Only needed for min and max temperature, as AgERA5 api provides daily accumulations for precipitation and reference ET

    # Prepare Copernicus Climate Data Store (CDS) API
    ensure_cds_dns()
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
                print("        *** DOWNLOADING CLIMATE DATA FROM AgERA5: " + variable + str(year) + " ***")
                retrieve_with_dns_fallback(
                    c,
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
            print("        *** DOWNLOADING SOIL MOISTURE DATA FROM ERA5-Land: " + variable + " ***")
            retrieve_with_dns_fallback(
                    c,
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