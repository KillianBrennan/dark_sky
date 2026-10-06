"""Regular LV95 analysis grids over the square domain (no drive-time mask)."""
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import geometry_mask

from . import config as C
from .terrain import load_dem


def isochrone_polys():
    iso = gpd.read_file(C.INTERIM / "isochrones.geojson").to_crs(C.CRS_M)
    return iso.sort_values("contour")


def make_grid(res):
    """Returns dict(transform, shape, inside mask, x, y, z (flattened over inside pixels),
    drive_min (15-min bins from the isochrone contours))."""
    from .fetch import domain_lv95
    iso = isochrone_polys()
    minx, miny, maxx, maxy = domain_lv95()
    w, h = int(round((maxx - minx) / res)), int(round((maxy - miny) / res))
    tr = rasterio.transform.from_origin(minx, maxy, res, res)
    inside = np.ones((h, w), bool)
    drive = np.full((h, w), np.nan, np.float32)   # informational: 15-min bins, NaN beyond 60 min
    for _, row in iso.sort_values("contour", ascending=False).iterrows():
        drive[~geometry_mask([row.geometry], (h, w), tr, all_touched=True)] = row.contour
    rows, cols = np.nonzero(inside)
    x = tr.c + (cols + 0.5) * res
    y = tr.f - (rows + 0.5) * res
    z, dtr = load_dem()
    fx = (x - dtr.c) / dtr.a - 0.5
    fy = (dtr.f - y) / -dtr.e - 0.5
    from scipy.ndimage import map_coordinates
    zz = map_coordinates(np.nan_to_num(z, nan=0.0), [fy, fx], order=1)
    return dict(transform=tr, shape=(h, w), inside=inside, rows=rows, cols=cols, x=x, y=y, z=zz,
                drive_min=drive[rows, cols], res=res)


def to_raster(g, values, fill=np.nan, dtype=np.float32):
    out = np.full(g["shape"], fill, dtype)
    out[g["rows"], g["cols"]] = values
    return out
