def soil(domain_path, res, basepath):

    from owslib.wcs import WebCoverageService
    import os
    import geopandas as gpd
    from soilgrids import SoilGrids
    import glob
    import xarray as xr
    import numpy as np
    import pdb # pdb.set_trace()
    from preproc_tools import makedirs


    ## Bounds from shapefile
    mask = gpd.read_file(domain_path)
    xmin, ymin, xmax, ymax = mask.total_bounds # in lat/lon
    xmin = float(np.floor(xmin * (1/res)) / (1/res)) # round bounds down to cell resolution
    ymin = float(np.floor(ymin * (1/res)) / (1/res)) # round bounds down to cell resolution
    xmax = float(np.ceil(xmax * (1/res)) / (1/res)) # round upper bounds up to cell resolution
    ymax = float(np.ceil(ymax * (1/res)) / (1/res)) # round upper bounds up to cell resolution

    w = round((xmax-xmin)/res)
    h = round((ymax-ymin)/res)


    ## Download soil data from ISRIC Soil grids
    target_dir = makedirs(basepath, 'rawdata', 'soilgrids')
    soil_grids = SoilGrids()

    # List of soilmaps to download
    maps =['clay', 'sand', 'silt', 'soc']

    print("        *** DOWNLOADING ISRIC SOILGRIDS DATA ***")

    # Browsing through input soil maps
    for i in maps:
        print(i)
        wcs = WebCoverageService('https://maps.isric.org/mapserv?map=/map/' + str(i) + '.map', version='2.0.1')

        # Choosing the mean layers
        names = [k for k in wcs.contents.keys() if k.endswith('mean')]

        # Download via soilgrids (https://github.com/gantian127/soilgrids)
        for layer in names:
            print(layer)
            outfile = os.path.join(target_dir, layer + '.tif')
            data = soil_grids.get_coverage_data(
                service_id= i,
                coverage_id= layer,
                crs= 'urn:ogc:def:crs:EPSG::4326',
                west= xmin,
                south= ymin,
                east= xmax,
                north= ymax,
                width= w,
                height= h,
                output=  outfile)


        ## Convert soil layers to mosaics for each depth layer
    print("Converting soil layers to .nc files")
    search_terms = ['*0-5*.tif', '*5-15*.tif', '*15-30*.tif', '*30-60*.tif', '*60-100*.tif', '*100-200*.tif'] # search files of different soil types
    
    for depth in search_terms:
        target_dir = makedirs(basepath, 'processed', '')
        targetfile = os.path.join(target_dir, 'soil_' + depth[1:-5] + '.nc')
    
        # Skip processing if file already exists
        if os.path.isfile(targetfile):
            print(f" Skipping {targetfile}, already exists.")
            continue
    
        file_to_mosaic = []
        path = os.path.join(basepath, 'rawdata', 'soilgrids', depth)
        ncs = glob.glob(path)
    
        for layer in ncs:
            src = xr.open_dataset(layer)
    
            # Rename bands as soil type, needs to be adjusted based on the data
            filename = os.path.basename(layer)
            soilname = filename.split('_')[0]
            src = src.rename({'band_data': soilname})
    
            file_to_mosaic.append(src)
    
        # Merge rasters
        mosaic = xr.merge(file_to_mosaic)
    

        # Convert clay, silt, and sand from (g/kg) to %.
        mosaic['Clay'] = mosaic.clay / 10
        mosaic['Sand'] = mosaic.sand / 10
        mosaic['Silt'] = mosaic.silt / 10
        #mosaic['Som'] = mosaic.soc * 1.724 / 100 # Soil organic matter derived from soil organic carbon (dg/kg) by multiplying with conventionally used factor of 1.724.
        mosaic['Som'] = mosaic.soc * 1.9 / 100 # EDIT: Pribyl 2010 shows that conventional factor of 1.724 is unfounded and too low and suggests 1.9 based on the median of values found in empirical studies.

        #mosaic = mosaic.drop(labels=['clay', 'soc', 'silt', 'sand'])
        mosaic = mosaic.drop_vars(['clay', 'soc', 'silt', 'sand'])

        # interpolating no data values
        som = mosaic.where(mosaic.Som != 0)
        mosaic = mosaic.where(mosaic.Clay != 0)
        mosaic['Clay'] = mosaic.Clay.rio.write_nodata(np.nan)
        mosaic['Clay'] = mosaic.Clay.rio.interpolate_na()
        mosaic['Silt'] = mosaic.Silt.rio.write_nodata(np.nan)
        mosaic['Silt'] = mosaic.Silt.rio.interpolate_na()
        mosaic['Sand'] = mosaic.Sand.rio.write_nodata(np.nan)
        mosaic['Sand'] = mosaic.Sand.rio.interpolate_na()
        som['Som'] = som.Som.rio.write_nodata(np.nan)
        mosaic['Som'] = som.Som.rio.interpolate_na()

        # Clipping with mask and save output files
        mosaic = mosaic.rio.clip(mask.geometry.values, mask.crs)
        mosaic = mosaic.drop_vars('band')
        #mosaic = mosaic.rename({'x': 'longitude', 'y': 'latitude'})

        target_dir = makedirs(basepath, 'processed', '')
        targetfile = os.path.join(target_dir, 'soil_' + depth[1:-5] + '.nc')
        mosaic.to_netcdf(targetfile)
