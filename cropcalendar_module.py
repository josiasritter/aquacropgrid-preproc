def cropcalendar(domain_path, basepath, referenceraster_path):

    import os
    import geopandas as gpd
    import glob
    import xarray as xr
    import rioxarray as rio
    from shapely.geometry import mapping
    from rasterio.warp import Resampling
    from preproc_tools import download_url, unzip_all, ensure_xy_dims, makedirs
    from xarray.coding.times import CFTimedeltaCoder

    # For reprojecting and clipping
    mask = gpd.read_file(domain_path)
    to_match = xr.open_dataset(referenceraster_path)
    to_match.rio.write_crs(4326, inplace=True)

    ## Download crop calendar dataset from GGCMI (Jägermeyr et al. 2021; https://zenodo.org/records/5062513)
    url = 'https://zenodo.org/api/records/5062513/files-archive'
    target_dir = makedirs(basepath, 'rawdata', 'cropcalendar')
    download_path = os.path.join(target_dir, 'ggcmi_cropcalendar.zip')
    unzipped_download_directory = download_path[:-4]  # assumes unzip goes to same-named folder

    #  Skip everything if already unzipped
    if os.path.exists(unzipped_download_directory):
        print(f" GGCMI crop calendar already unzipped, skipping: {unzipped_download_directory}")
    else:
        #  Skip download if ZIP exists
        if os.path.exists(download_path):
            print(f" ZIP already exists, skipping download: {download_path}")
        else:
            print(" Downloading GGCMI crop calendar data")
            print('   URL:', url)
            download_url(url, download_path=download_path)

        print("Unzipping GGCMI crop calendar data...")
        unzip_all(dir=target_dir)


    ## Convert separate files to rasterstack, reproject (nearest neighbour) to 0.05 degrees, and clip to domain

    # Dictionary that connects GGCMI crop ID's to crop names in AquaCrop. Used for layer naming
    crop_dict = {'bar': 'Barley', 'cot': 'Cotton', 'bea': 'DryBean', 'mai': 'Maize', 'ri1': 'PaddyRice1', 'ri2': 'PaddyRice2', 'pot': 'Potato', 'sor': 'Sorghum', 'soy': 'Soybean', 'sgb': 'SugarBeet', 'sgc': 'SugarCane', 'sun': 'Sunflower', 'swh': 'Wheat_summer', 'wwh': 'Wheat_winter', 'cas': 'Cassava'}

    search_criteria = os.path.join(unzipped_download_directory, '*.nc4')
    layers = glob.glob(search_criteria)
    file_to_mosaic = []
    for layer in layers:

        # Get crop ID from layer name and preserve only if crop type is supported by AquaCrop
        cropID = layer[-43:-40]
        technique = layer[-39:-37]
        if cropID not in crop_dict:
            os.remove(layer)
            continue

        src = xr.open_dataset(layer, decode_timedelta=False)       
        src = ensure_xy_dims(src)
        #src = src.drop_vars(['growing_season_length', 'data_source_used'])
        src = src.drop_vars(['maturity_day', 'data_source_used'])
        src.rio.write_crs(4326, inplace=True)

        # Reproject and clip to model domain
        src_proj = src.rio.reproject_match(to_match)#, resampling=Resampling.nearest)  # nearest-neighbor
        src_clip = src_proj.rio.clip(mask.geometry.apply(mapping))

        # Rename variables from FAO crop ID to crop strings used by AquaCrop
        layername_base = crop_dict.get(cropID) + '_' + technique
        src_clip[layername_base + '_planting'] = src_clip['planting_day']
        #src_clip[layername_base + '_harvest'] = src_clip['maturity_day']
        src_clip[layername_base + '_growing_season_length'] = src_clip['growing_season_length']

        # Select relevant data variables for merging
        file_to_mosaic.append(src_clip[layername_base + '_planting'])
        #file_to_mosaic.append(src_clip[layername_base + '_harvest'])
        file_to_mosaic.append(src_clip[layername_base + '_growing_season_length'])

    # Merge data into one file
    src_mosaic = xr.merge(file_to_mosaic)
    target_dir = makedirs(basepath, 'processed', '')
    targetfile = os.path.join(target_dir, 'cropcalendar.nc')
    src_mosaic.to_netcdf(targetfile)
