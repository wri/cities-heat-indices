from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, timedelta
from typing import Optional, Tuple, List

import ee

from city_metrix.metrix_model import Layer, GeoExtent


GTIFF_FILE_EXTENSION = ".tif"
DEFAULT_SPATIAL_RESOLUTION = 30  # Landsat native
DEFAULT_HALF_WINDOW_DAYS = 45    # 90-day window


@dataclass
class LstSeasonWindow:
    year_offset: int
    start: str  # YYYY-MM-DD
    end: str    # YYYY-MM-DD


class LandsatLstP95(Layer):
    """
    Pixel-wise 95th percentile LST (°C) across a 3-window seasonal stack.

    Window logic (each window is ~3 months centered on center_date):
      - Normally: year-1, year, year+1
      - If year+1 window is not complete as of today (UTC): year-2, year-1, year

    Output band: LST_C_p95
    """
    OUTPUT_FILE_FORMAT = GTIFF_FILE_EXTENSION
    MAJOR_NAMING_ATTS = ["center_date", "stat", "half_window_days"]
    MINOR_NAMING_ATTS = None
    PROCESSING_TILE_SIDE_M = 5000  # similar to ESA; adjust if you want
    STAT = "p95"

    def __init__(
        self,
        center_date: date,
        percentile: int = 95,
        half_window_days: int = DEFAULT_HALF_WINDOW_DAYS,
        include_l9: bool = True,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.center_date = center_date
        self.percentile = percentile
        self.half_window_days = half_window_days
        self.include_l9 = include_l9

    # ---------- internal helpers ----------
    def _centered_window(self, year_offset: int) -> Tuple[str, str]:
        c = date(self.center_date.year + year_offset, self.center_date.month, self.center_date.day)
        start = c - timedelta(days=self.half_window_days)
        end = c + timedelta(days=self.half_window_days)
        return start.isoformat(), end.isoformat()

    def _choose_offsets(self) -> List[int]:
        today = datetime.utcnow().date()
        _, end_plus1 = self._centered_window(+1)
        end_plus1_dt = datetime.strptime(end_plus1, "%Y-%m-%d").date()
        return [-1, 0, 1] if end_plus1_dt <= today else [-2, -1, 0]

    @staticmethod
    def _cloud_mask_l2(img: ee.Image) -> ee.Image:
        # QA_PIXEL bits: 3=cloud, 4=cloud shadow
        qa = img.select("QA_PIXEL")
        cloud = qa.bitwiseAnd(1 << 3).neq(0)
        shadow = qa.bitwiseAnd(1 << 4).neq(0)
        return img.updateMask(cloud.Or(shadow).Not())

    @staticmethod
    def _st_b10_to_celsius(img: ee.Image) -> ee.Image:
        # Landsat C2 L2 ST_B10 scale/offset -> Kelvin -> Celsius
        lst_c = (
            img.select("ST_B10")
            .multiply(0.00341802)
            .add(149.0)
            .subtract(273.15)
            .rename("LST_C")
        )
        return img.addBands(lst_c, overwrite=True)

    def _build_lst_collection(self, aoi_geom: ee.Geometry, start: str, end: str) -> ee.ImageCollection:
        l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
        ic = l8
        if self.include_l9:
            l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
            ic = ic.merge(l9)

        return (
            ic.filterBounds(aoi_geom)
              .filterDate(start, end)
              .select(["ST_B10", "QA_PIXEL"])
              .map(self._cloud_mask_l2)
              .map(self._st_b10_to_celsius)
        )

    # ---------- public API ----------
    def get_windows(self) -> List[LstSeasonWindow]:
        offsets = self._choose_offsets()
        out: List[LstSeasonWindow] = []
        for off in offsets:
            s, e = self._centered_window(off)
            out.append(LstSeasonWindow(off, s, e))
        return out

    def get_data(
        self,
        bbox: GeoExtent,
        spatial_resolution: int = DEFAULT_SPATIAL_RESOLUTION,
        resampling_method=None,
    ) -> ee.Image:
        if resampling_method is not None:
            raise Exception("resampling_method can not be specified.")
        spatial_resolution = DEFAULT_SPATIAL_RESOLUTION if spatial_resolution is None else spatial_resolution

        # ee_rectangle = bbox.to_ee_rectangle()
        # aoi_geom = ee_rectangle  # geometry works fine here
        minx, miny, maxx, maxy = bbox.bbox if hasattr(bbox, "bbox") else bbox.bounds
        aoi_geom = ee.Geometry.Rectangle([float(minx), float(miny), float(maxx), float(maxy)],
                                         proj="EPSG:4326", geodesic=False)

        merged = ee.ImageCollection([])
        for w in self.get_windows():
            merged = merged.merge(self._build_lst_collection(aoi_geom, w.start, w.end))

        # Pixel-wise percentile across time
        p_img = (
            merged.select("LST_C")
                  .reduce(ee.Reducer.percentile([self.percentile]))
                  .rename(f"LST_C_p{self.percentile}")
                  .clip(aoi_geom)
        )

        # Keep projection reasonable (native Landsat)
        # Setting a default projection here helps downloads be consistent.
        proj = merged.first().select("LST_C").projection()
        p_img = p_img.setDefaultProjection(proj)

        return p_img


# #!/usr/bin/env python3
# import argparse
# import os
# from datetime import datetime, timedelta, date
# 
# import ee
# import geemap  # pip install geemap
# 
# import urllib.request
# import json
# 
# HALF_WINDOW_DAYS = 45  # 90-day window centered on date
# 
# def parse_args():
#     p = argparse.ArgumentParser()
#     p.add_argument("--date", required=True, help="Center date YYYY-MM-DD")
#     p.add_argument("--boundary", required=True, help="Path to GeoJSON OR EE FeatureCollection asset id")
#     p.add_argument("--out", required=True, help="Output GeoTIFF path, e.g. /tmp/lst_p95.tif")
#     p.add_argument("--scale", type=int, default=30)
#     p.add_argument("--max_pixels", type=float, default=1e13)
#     return p.parse_args()
# 
# def centered_window(center: date, year_offset: int) -> tuple[str, str]:
#     c = date(center.year + year_offset, center.month, center.day)
#     start = c - timedelta(days=HALF_WINDOW_DAYS)
#     end = c + timedelta(days=HALF_WINDOW_DAYS)
#     return start.isoformat(), end.isoformat()
# 
# def choose_offsets(center: date) -> list[int]:
#     # Normally [-1,0,+1]; if +1 window not yet complete, use [-2,-1,0]
#     today = datetime.utcnow().date()
#     _, end_plus1 = centered_window(center, +1)
#     end_plus1_dt = datetime.strptime(end_plus1, "%Y-%m-%d").date()
#     return [-1, 0, 1] if end_plus1_dt <= today else [-2, -1, 0]
# 
# def load_aoi(boundary: str) -> ee.Geometry:
#     # URL GeoJSON
#     if boundary.startswith("http://") or boundary.startswith("https://"):
#         with urllib.request.urlopen(boundary) as resp:
#             gj = json.loads(resp.read().decode("utf-8"))
# 
#         if gj.get("type") == "FeatureCollection":
#             return ee.FeatureCollection(gj).geometry()
#         else:
#             return ee.Geometry(gj)
# 
#     # Local file GeoJSON
#     if os.path.exists(boundary):
#         with open(boundary, "r", encoding="utf-8") as f:
#             gj = json.load(f)
#         if gj.get("type") == "FeatureCollection":
#             return ee.FeatureCollection(gj).geometry()
#         else:
#             return ee.Geometry(gj)
# 
#     # EE asset id (FeatureCollection)
#     return ee.FeatureCollection(boundary).geometry()
# 
# def cloud_mask_l2(img: ee.Image) -> ee.Image:
#     qa = img.select("QA_PIXEL")
#     cloud  = qa.bitwiseAnd(1 << 3).neq(0)
#     shadow = qa.bitwiseAnd(1 << 4).neq(0)
#     return img.updateMask(cloud.Or(shadow).Not())
# 
# def st_b10_to_celsius(img: ee.Image) -> ee.Image:
#     lst_c = (img.select("ST_B10")
#              .multiply(0.00341802)
#              .add(149.0)
#              .subtract(273.15)
#              .rename("LST_C"))
#     return img.addBands(lst_c, overwrite=True)
# 
# def build_lst_collection(aoi: ee.Geometry, start: str, end: str) -> ee.ImageCollection:
#     l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
#     # l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
#     return (l8.filterBounds(aoi)
#             .filterDate(start, end)
#             .select(["ST_B10", "QA_PIXEL"])
#             .map(cloud_mask_l2)
#             .map(st_b10_to_celsius))
# 
# def main():
#     args = parse_args()
# 
#     # First run (one time on your machine): ee.Authenticate()
#     ee.Initialize()
# 
#     center = datetime.strptime(args.date, "%Y-%m-%d").date()
#     aoi = load_aoi(args.boundary)
# 
#     offsets = choose_offsets(center)
# 
#     merged = ee.ImageCollection([])
# 
#     windows = []
#     for off in offsets:
#         s, e = centered_window(center, off)
#         windows.append((off, s, e))
#         merged = merged.merge(build_lst_collection(aoi, s, e))
# 
#     # Pixel-wise 95th percentile across time
#     p95_img = (merged
#                .select("LST_C")
#                .reduce(ee.Reducer.percentile([95]))
#                .rename("LST_C_p95")
#                .clip(aoi))
# 
#     # (Optional) print a simple AOI summary of the p95 image
#     # stats = p95_img.reduceRegion(
#     #     reducer=ee.Reducer.median().combine(ee.Reducer.mean(), sharedInputs=True),
#     #     geometry=aoi,
#     #     scale=args.scale,
#     #     maxPixels=args.max_pixels,
#     #     bestEffort=True,
#     # ).getInfo()
# 
#     print("Input date:", args.date)
#     print("Used 90-day windows (±45 days):")
#     for off, s, e in windows:
#         print(f"  offset {off:+}: {s} → {e}")
#     # print("AOI summary of pixel-wise p95 image (°C):", stats)
# 
#     # Download GeoTIFF locally (tiles large images when needed)
#     os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
#     geemap.download_ee_image(
#         image=p95_img,
#         filename=args.out,
#         region=aoi
#         # scale=args.scale
#         # crs=None,          # let EE choose unless you want a specific CRS
#         # max_pixels=args.max_pixels
#     )
# 
#     print("Wrote:", args.out)
# 
# if __name__ == "__main__":
#     main()
