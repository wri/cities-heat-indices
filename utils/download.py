import os
import urllib.request
from datetime import datetime

import ee
from city_metrix.metrix_model import GeoExtent

from .LST import *
from city_metrix.layers import FractionalVegetationPercent
from utils.seasonal_albedo import SeasonalAlbedo


def download_ee_image_geotiff(
    image: ee.Image,
    out_tif: str,
    region: ee.Geometry,
    scale: int = 30,
    crs: str | None = None,
):
    """
    Download a GeoTIFF via EE's getDownloadURL (good for tile-sized regions).
    """
    params = {
        "scale": scale,
        "region": region,
        "format": "GEO_TIFF",
        "filePerBand": False,
    }
    if crs is not None:
        params["crs"] = crs

    url = image.getDownloadURL(params)
    os.makedirs(os.path.dirname(os.path.abspath(out_tif)), exist_ok=True)
    urllib.request.urlretrieve(url, out_tif)


def get_lst(city, bbox: GeoExtent, grid_cell_id, data_path, center_date: str, copy_to_s3=False, crs=None):
    """
    Fetch Landsat LST pixel-wise p95 for one grid cell and write a GeoTIFF.

    Args:
        city (str)
        bbox (GeoExtent)
        grid_cell_id (int|str)
        data_path (str)
        center_date (str): "YYYY-MM-DD"
        copy_to_s3 (bool)
        crs (str|None): optional EPSG like "EPSG:4326" or an EE CRS string
    """
    ee.Initialize()

    lst_path = f"{data_path}/{city}/lst"
    out_file = f"{lst_path}/lst_p95_{center_date}_{grid_cell_id}.tif"

    if os.path.exists(out_file):
        print(f"LST already exists at {out_file}, skipping.")
        return out_file

    print(f"Fetching LST p95 for {city}, cell {grid_cell_id}...")

    center_dt = datetime.strptime(center_date, "%Y-%m-%d").date()
    layer = LandsatLstP95(center_date=center_dt, percentile=95, half_window_days=45, include_l9=True)

    img = layer.get_data(bbox, spatial_resolution=30)  # returns ee.Image band "LST_C_p95"
    # region = bbox.to_ee_rectangle()
    if hasattr(bbox, "bbox"):
        minx, miny, maxx, maxy = bbox.bbox
    else:
        minx, miny, maxx, maxy = bbox.bounds
    minx, miny, maxx, maxy = map(float, [minx, miny, maxx, maxy])

    # Build a valid EE rectangle in WGS84
    region = ee.Geometry.Rectangle(
        [minx, miny, maxx, maxy],
        proj="EPSG:4326",
        geodesic=False
    )

    os.makedirs(lst_path, exist_ok=True)
    download_ee_image_geotiff(
        image=img,
        out_tif=out_file,
        region=region,
        scale=30,
        crs=crs,
    )

    print("Wrote:", out_file)

    if copy_to_s3:
        to_s3(out_file, data_path)  # assuming you already have this utility

    return out_file

def get_fr(city, bbox: GeoExtent, year, grid_cell_id, data_path, copy_to_s3=False):
    """
    Fetch fractional vegetation (Fr) for one grid cell and write a GeoTIFF.

    Args:
        city (str)
        bbox (GeoExtent)
        year (int): year of imagery (e.g. 2024)
        grid_cell_id (int|str)
        data_path (str)
        copy_to_s3 (bool)
    """
    ee.Initialize()
    
    fr_path = f"{data_path}/{city}/fr"
    fr_file = f"{fr_path}/fr_{year}_{grid_cell_id}.tif"

    # skip if file exists
    if os.path.exists(fr_file):
        print(f"Fr data already exists at {fr_file}, skipping fetch.")
        return fr_file

    print(f"Fetching Fr data for {city}, cell {grid_cell_id}…")
    ee.Initialize()

    # Create the FractionalVegetationPercent layer for the given year
    fr_layer = FractionalVegetationPercent(year=year)
    # get_data returns an xarray DataArray with band "Fr"
    fr_data = fr_layer.get_data(bbox)

    # ensure output directory exists
    os.makedirs(fr_path, exist_ok=True)
    # write GeoTIFF
    fr_data.rio.to_raster(fr_file)
    print("Wrote:", fr_file)

    if copy_to_s3:
        to_s3(fr_file, data_path)

    return fr_file

from utils.seasonal_albedo import SeasonalAlbedo

def get_albedo(city, bbox, grid_cell_id, data_path, year=None, copy_to_s3=False):
    """Fetch seasonal albedo for one grid cell and write a GeoTIFF."""
    ee.Initialize()

    alb_path = f"{data_path}/{city}/albedo"
    out_file = f"{alb_path}/albedo_{year or 'recent'}_{grid_cell_id}.tif"

    if os.path.exists(out_file):
        print(f"Albedo already exists at {out_file}, skipping.")
        return out_file

    print(f"Fetching albedo for {city}, cell {grid_cell_id}…")

    layer = SeasonalAlbedo(year=year)
    img = layer.get_data(bbox, spatial_resolution=10)

    # build region geometry from bbox
    if hasattr(bbox, "bbox"):
        minx, miny, maxx, maxy = bbox.bbox
    else:
        minx, miny, maxx, maxy = bbox.bounds
    region = ee.Geometry.Rectangle([float(minx), float(miny), float(maxx), float(maxy)],
                                   proj="EPSG:4326", geodesic=False)

    os.makedirs(alb_path, exist_ok=True)
    download_ee_image_geotiff(
        image=img.select("albedo"),
        out_tif=out_file,
        region=region,
        scale=10,
        crs=None,
    )

    print("Wrote:", out_file)
    if copy_to_s3:
        to_s3(out_file, data_path)

    return out_file

