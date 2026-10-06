"""VIIRS DNB source term (EOG VNL v2.2 annual, average_masked, nW/cm^2/sr).

EOG requires a (free) login, so the file is not fetched automatically: download it from
https://eogdata.mines.edu/nighttime_light/annual/v22/<year>/ into data/raw/viirs/ and
set config.VIIRS_FILE. Any GeoTIFF of DNB radiance in nW/cm^2/sr covering the domain
+ SOURCE_RADIUS_KM works (e.g. a GEE export of NOAA/VIIRS/DNB/ANNUAL_V22).

Radiance is treated as an emission *source term* proportional to the nadir upward
intensity of each pixel (radiance x pixel area), not as sky brightness.
"""
import math

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds

from . import config as C
from .fetch import domain_bounds

NEAR_RES = 1000.0
FAR_RES = 5000.0


def clip():
    out = C.INTERIM / "viirs_clip.tif"
    if out.exists():
        return out
    minx, miny, maxx, maxy = domain_bounds()
    lat0 = (miny + maxy) / 2
    dlat = C.SOURCE_RADIUS_KM / 111.0 + 0.1
    dlon = C.SOURCE_RADIUS_KM / (111.0 * math.cos(math.radians(lat0 + dlat))) + 0.1
    src_path = C.VIIRS_FILE
    path = f"/vsigzip/{src_path}" if str(src_path).endswith(".gz") else str(src_path)
    with rasterio.open(path) as src:
        win = from_bounds(minx - dlon, miny - dlat, maxx + dlon, maxy + dlat, src.transform).round_offsets().round_lengths()
        arr = src.read(1, window=win)
        tr = src.window_transform(win)
        prof = dict(driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1, dtype="float32",
                    crs=src.crs, transform=tr, compress="deflate", tiled=True)
    with rasterio.open(out, "w", **prof) as dst:
        dst.write(arr.astype(np.float32), 1)
    return out


def sources(led_factor=1.0):
    """Returns (near, far) source sets on LV95 grids, each a dict of x, y, w (relative
    nadir intensity, nW/sr * 1e-9 scaled), and for `far` the index ranges into `near`
    of the 1 km cells each 5 km block contains (so close blocks can be resolved finely)."""
    with rasterio.open(clip()) as src:
        rad = src.read(1).astype(np.float64)
        tr = src.transform
    rad = np.where(np.isfinite(rad) & (rad > C.VIIRS_MIN_RADIANCE), rad, 0.0)
    rows, cols = np.nonzero(rad)
    lon = tr.c + (cols + 0.5) * tr.a
    lat = tr.f + (rows + 0.5) * tr.e
    # pixel area (cm^2) on the sphere
    area = (abs(tr.a) * math.pi / 180 * C.R_EARTH) * (abs(tr.e) * math.pi / 180 * C.R_EARTH) * np.cos(np.radians(lat)) * 1e4
    w = rad[rows, cols] * area * 1e-9 * 1e-9 * led_factor   # nW/sr -> W/sr, then 1e-9 for O(1) numbers
    x, y = Transformer.from_crs(C.CRS_LL, C.CRS_M, always_xy=True).transform(lon, lat)

    # elevation of each source from the DEM where available (lowland default elsewhere)
    from .terrain import load_dem
    z, dtr = load_dem()
    c = ((x - dtr.c) / dtr.a).astype(int)
    r = ((y - dtr.f) / dtr.e).astype(int)
    ok = (r >= 0) & (r < z.shape[0]) & (c >= 0) & (c < z.shape[1])
    zs = np.full(x.shape, C.SOURCE_ALT)
    zs[ok] = np.nan_to_num(z[r[ok], c[ok]], nan=C.SOURCE_ALT)

    def agg(res, x, y, w, zs):
        kx = np.floor(x / res).astype(np.int64)
        ky = np.floor(y / res).astype(np.int64)
        key = kx * 10_000_000 + ky
        uk, inv = np.unique(key, return_inverse=True)
        W = np.bincount(inv, weights=w)
        X = np.bincount(inv, weights=w * x) / W            # intensity-weighted centroid
        Y = np.bincount(inv, weights=w * y) / W
        Z = np.bincount(inv, weights=w * zs) / W
        return uk, inv, X, Y, Z, W

    uk1, _, X1, Y1, Z1, W1 = agg(NEAR_RES, x, y, w, zs)
    # sort 1 km cells by their parent 5 km block so each block maps to a contiguous range
    kx1 = uk1 // 10_000_000
    ky1 = uk1 % 10_000_000
    ratio = int(FAR_RES / NEAR_RES)
    parent = (np.floor_divide(kx1, ratio)) * 10_000_000 + np.floor_divide(ky1, ratio)
    order = np.argsort(parent, kind="stable")
    X1, Y1, Z1, W1, parent = X1[order], Y1[order], Z1[order], W1[order], parent[order]
    up, start, counts = np.unique(parent, return_index=True, return_counts=True)
    W5 = np.add.reduceat(W1, start)
    X5 = np.add.reduceat(W1 * X1, start) / W5
    Y5 = np.add.reduceat(W1 * Y1, start) / W5
    Z5 = np.add.reduceat(W1 * Z1, start) / W5
    near = dict(x=X1, y=Y1, z=Z1, w=W1)
    far = dict(x=X5, y=Y5, z=Z5, w=W5, start=start.astype(np.int64), count=counts.astype(np.int64))
    return near, far
