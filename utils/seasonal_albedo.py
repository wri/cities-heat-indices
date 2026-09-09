"""
Seasonal albedo computation using Sentinel‑2.

This module provides a helper to compute a seasonal median albedo composite
from Sentinel‑2 surface reflectance imagery. By default, it computes the
most recent complete growing‑season composite: for the Northern Hemisphere
this is the June–August (JJA) season and for the Southern Hemisphere the
December–February (DJF) season. If you supply a specific ``year``
argument, the composite will instead span the JJA season of that year in
the north or the DJF season labelled by that year in the south (i.e.,
December of ``year-1`` through March of ``year``).

Example usage::

    import ee
    from seasonal_albedo import seasonal_s2_albedo

    ee.Initialize()
    geom = ee.Geometry.Point([-122.335167, 47.608013]).buffer(10000)
    # Most recent complete season
    img = seasonal_s2_albedo(geom)
    # JJA 2024 in the north (or DJF 2023/2024 in the south)
    img_2024 = seasonal_s2_albedo(geom, year=2024)

The returned :class:`ee.Image` has one band named ``"albedo"`` and
properties ``"season_start"``, ``"season_end"`` and ``"hemisphere"``.
"""

import ee
from datetime import datetime, date

# Cloud Score+ threshold used to mask clouds and haze. Values in the
# cloud score band range from 0–1; pixels with score below the
# threshold are masked out. Thresholds between 0.50 and 0.65 work well.
CLEAR_THRESHOLD = 0.60


def _add_albedo(image: ee.Image) -> ee.Image:
    """Compute the Bonafoni albedo for a Sentinel‑2 image.

    Parameters
    ----------
    image : ee.Image
        A Sentinel‑2 surface reflectance image in the
        ``COPERNICUS/S2_SR_HARMONIZED`` collection. Must have the
        reflectance bands B2 (blue), B3 (green), B4 (red), B8 (NIR),
        B11 (SWIR1) and B12 (SWIR2).

    Returns
    -------
    ee.Image
        The input image with an additional ``"albedo"`` band computed
        using the linear coefficients of Bonafoni et al. (2020).
    """
    image = ee.Image(image)
    albedo = image.expression(
        "B*Bw + G*Gw + R*Rw + NIR*NIRw + SW1*SW1w + SW2*SW2w",
        {
            "B": image.select("B2"),
            "G": image.select("B3"),
            "R": image.select("B4"),
            "NIR": image.select("B8"),
            "SW1": image.select("B11"),
            "SW2": image.select("B12"),
            "Bw": 0.2266,
            "Gw": 0.1236,
            "Rw": 0.1573,
            "NIRw": 0.3417,
            "SW1w": 0.1170,
            "SW2w": 0.0338,
        },
    ).rename("albedo")
    return image.addBands(albedo)


def _season_window_for(geom: ee.Geometry, year: int | None = None) -> dict:
    """Compute the start and end of the relevant seasonal window.

    For the Northern Hemisphere the growing season is June–August (JJA);
    for the Southern Hemisphere it is December–February (DJF). If
    ``year`` is provided, the window will correspond to that year: for
    JJA the window starts on 1 June ``year`` and ends on 1 September
    ``year``; for DJF the window starts on 1 December ``year-1`` and
    ends on 1 March ``year``. If ``year`` is ``None``, the most recent
    complete season relative to today’s date is returned.

    Parameters
    ----------
    geom : ee.Geometry
        The geometry of interest used to determine hemisphere.
    year : int | None, optional
        The calendar year of the season; see above. If ``None``, the
        function falls back to the most recent complete season.

    Returns
    -------
    dict
        A dictionary containing the keys ``"start"`` and ``"end"`` with
        :class:`ee.Date` values, and ``"hemisphere"`` with a string
        (``"north"`` or ``"south"``).
    """
    g = ee.Geometry(geom)
    # Determine hemisphere from centroid latitude
    centroid = g.centroid(1).coordinates().getInfo()
    lat = centroid[1]
    is_north = lat >= 0

    if year is None:
        # Determine the year based on today and whether the current
        # seasonal window has completed
        today = datetime.utcnow().date()
        yr = today.year
        if is_north:
            # June 1 to Sep 1 of the most recent complete year
            start_this = date(yr, 6, 1)
            end_this = date(yr, 9, 1)
            if today >= end_this:
                start = start_this
                end = end_this
            else:
                start = date(yr - 1, 6, 1)
                end = date(yr - 1, 9, 1)
        else:
            # Dec 1 to Mar 1 of the most recent complete season
            mar1_this = date(yr, 3, 1)
            if today >= mar1_this:
                # last complete DJF spans Dec of previous year to Mar this year
                start = date(yr - 1, 12, 1)
                end = date(yr, 3, 1)
            else:
                # last complete DJF spans Dec of two years ago to Mar of last year
                start = date(yr - 2, 12, 1)
                end = date(yr - 1, 3, 1)
    else:
        # Use the provided year to define the season
        if is_north:
            # JJA: Jun–Aug of the specified year
            start = date(year, 6, 1)
            end = date(year, 9, 1)
        else:
            # DJF: Dec of year-1 to Mar of year
            start = date(year - 1, 12, 1)
            end = date(year, 3, 1)

    hemi = "north" if is_north else "south"
    start_ee = ee.Date(start.isoformat())
    end_ee = ee.Date(end.isoformat())
    return {"start": start_ee, "end": end_ee, "hemisphere": hemi}


def seasonal_s2_albedo(geom: ee.Geometry | ee.Feature | ee.FeatureCollection, year: int | None = None) -> ee.Image:
    """Create a Sentinel‑2 seasonal albedo composite for a geometry.

    This function returns a median composite of albedo over the
    appropriate growing season defined by the geometry’s hemisphere and
    an optional year. If ``year`` is ``None``, it uses the most recent
    complete season. Otherwise, it uses the JJA season of ``year`` in
    the north or the DJF season labelled by ``year`` in the south.

    Parameters
    ----------
    geom : ee.Geometry or ee.Feature or ee.FeatureCollection
        The area of interest. The albedo image will be clipped to this
        geometry’s bounds.
    year : int | None, optional
        The target year for the season. See :func:`_season_window_for`.

    Returns
    -------
    ee.Image
        An image with a single band ``"albedo"`` clipped to the input
        geometry. The image has properties ``"season_start"`` and
        ``"season_end"`` formatted as ``YYYY-MM-dd``, and
        ``"hemisphere"``.
    """
    # Normalize input to an EE geometry
    if isinstance(geom, ee.feature.Feature):
        g = geom.geometry()
    elif isinstance(geom, ee.featurecollection.FeatureCollection):
        g = geom.geometry()
    else:
        g = ee.Geometry(geom)

    # Determine the seasonal window and hemisphere
    win = _season_window_for(g, year)
    start = win["start"]
    end = win["end"]
    hemi = win["hemisphere"]

    # Sentinel‑2 surface reflectance collection and cloud score+ dataset
    S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    S2CS = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")

    # Filter images by location and season, join with cloud score+, mask
    # clouds/haze using the CLEAR_THRESHOLD, and scale reflectance to [0,1]
    s2_masked = (
        S2.filterBounds(g)
        .filterDate(start, end)
        .linkCollection(S2CS, ["cs"])
        .map(lambda img: ee.Image(img)
             .updateMask(ee.Image(img).select("cs").gte(CLEAR_THRESHOLD))
             .divide(10000))
    )

    # Add albedo band to each image and select it
    albedo_ic = s2_masked.map(_add_albedo).select("albedo")

    # Compute a median composite; if no images are available, return a
    # masked zero image as a placeholder
    alb = ee.Image(
        ee.Algorithms.If(
            albedo_ic.size().gt(0),
            albedo_ic.median(),
            ee.Image(0).rename("albedo").updateMask(ee.Image(0)),
        )
    )

    # Clip to the geometry bounds with a small buffer (1 m)
    alb = alb.clip(g.bounds(1))

    # Ensure values are within [0,1]
    alb = alb.where(alb.gte(1), 1).where(alb.lte(0), 0)

    # Attach metadata
    alb = alb.set({
        "season_start": start.format("YYYY-MM-dd"),
        "season_end": end.format("YYYY-MM-dd"),
        "hemisphere": hemi,
    })

    return alb
  
import ee
from datetime import datetime, date

# existing helpers: CLEAR_THRESHOLD, _add_albedo, _season_window_for, seasonal_s2_albedo …

try:
    from city_metrix.metrix_model import Layer, GeoExtent
except Exception:
    # Fallback so the module can be imported outside of city_metrix
    Layer = object  # type: ignore
    GeoExtent = object  # type: ignore

class SeasonalAlbedo(Layer):
    """Compute a seasonal albedo composite for a bounding box.

    Parameters
    ----------
    year : int | None, optional
        Target year for the season.  If None (default), uses the most
        recent complete season: JJA in the north or DJF in the south.
    clear_threshold : float, optional
        Cloud Score+ threshold (0–1).  Lower values mask more clouds.
    **kwargs : dict
        Passed to the parent Layer constructor.
    """
    OUTPUT_FILE_FORMAT = ".tif"
    MAJOR_NAMING_ATTS = ["year"]
    MINOR_NAMING_ATTS = None
    PROCESSING_TILE_SIDE_M = 5000

    def __init__(self, year: int | None = None, clear_threshold: float = CLEAR_THRESHOLD, **kwargs):
        super().__init__(**kwargs)
        self.year = year
        self.clear_threshold = clear_threshold

    def get_data(self,
                 bbox: GeoExtent,
                 spatial_resolution: int = 10,
                 resampling_method: str | None = None) -> ee.Image:
        """Return a seasonal albedo image clipped to the bounding box.

        bbox : GeoExtent
            A lon/lat bounding box; GeoExtent may have a `.bbox` or `.bounds`
            attribute.  If those aren’t present, it’s treated as an iterable
            (minx, miny, maxx, maxy).
        spatial_resolution : int, optional
            Pixel size in metres.  Sentinel‑2 data is native at 10 m.
        resampling_method : str | None, optional
            Not used; specifying it raises an exception.
        """
        if resampling_method is not None:
            raise Exception("resampling_method can not be specified.")

        # extract numeric bounds
        try:
            if hasattr(bbox, "bbox"):
                minx, miny, maxx, maxy = bbox.bbox
            else:
                minx, miny, maxx, maxy = bbox.bounds
        except Exception:
            minx, miny, maxx, maxy = bbox  # fallback

        minx, miny, maxx, maxy = map(float, [minx, miny, maxx, maxy])
        aoi = ee.Geometry.Rectangle([minx, miny, maxx, maxy], proj="EPSG:4326", geodesic=False)

        # compute albedo; self.year may be None or an int
        img = seasonal_s2_albedo(aoi, year=self.year)

        # set default projection for consistent downloads
        try:
            proj = img.select("albedo").projection()
            img = img.setDefaultProjection(proj)
        except Exception:
            pass

        return img.clip(aoi)

