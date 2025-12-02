import os
import geopandas as gpd
from datetime import datetime

def validate_inputs(domain_path, start_year, end_year, api_token):
    # --- 1. Check that the input file exists ---
    if not os.path.exists(domain_path):
        raise FileNotFoundError(f"❌ The input file '{domain_path}' does not exist.")

    # --- 2. Try reading the file with GeoPandas ---
    try:
        gdf = gpd.read_file(domain_path)
    except Exception as e:
        raise ValueError(f"❌ Could not read file '{domain_path}'. Error: {e}")

    # --- 3. Check that the file contains polygon geometries ---
    geom_types = gdf.geom_type.unique()
    allowed_types = {"Polygon", "MultiPolygon"}

    if not all(g in allowed_types for g in geom_types):
        raise ValueError(f"❌ The file must contain only Polygon or MultiPolygon geometries. Found: {geom_types}")

    # --- 4. Check CRS (Coordinate Reference System) ---
    if gdf.crs is None:
        raise ValueError("❌ The file has no CRS defined. It must be WGS84 (EPSG:4326).")

    epsg_code = gdf.crs.to_epsg()
    if epsg_code != 4326:
        raise ValueError(f"❌ The CRS of the file is {gdf.crs}, not WGS84 (EPSG:4326).")

    # --- 5. Print the extents (in lat/lon) ---
    minx, miny, maxx, maxy = gdf.total_bounds
    print(f"✅ Polygon extents (lat/lon):")
    print(f"   Longitude: {minx:.4f} to {maxx:.4f}")
    print(f"   Latitude:  {miny:.4f} to {maxy:.4f}")

    # --- Additional warning for small AOI ---
    dx = abs(maxx - minx)
    dy = abs(maxy - miny)

    if dx < 0.25 or dy < 0.25:
        print("⚠️  Warning: The spatial extent of the area of interest is very small.")
        print("   Climate datasets have a coarse resolution (0.25° in the case of ReferenceET).")
        print("   The code may fail or produce unreliable results due to low spatial resolution of climate data.")

    # --- 6. Check year values ---
    current_year = datetime.now().year
    for name, value in [("start_year", start_year), ("end_year", end_year)]:
        if not isinstance(value, (int, float)):
            raise TypeError(f"❌ {name} must be a numeric value.")
        if not (1950 <= value <= current_year):
            raise ValueError(f"❌ {name} must be between 1950 and {current_year}. Given: {value}")

    if start_year > end_year:
        raise ValueError(f"❌ start_year ({start_year}) cannot be greater than end_year ({end_year}).")

    # --- 7. Check API token ---
    if not isinstance(api_token, str):
        raise TypeError("❌ API token must be a string.")
    if len(api_token.strip()) <= 30:
        raise ValueError("❌ API token must be longer than 30 characters.")

    print("✅ All input checks passed successfully.")
