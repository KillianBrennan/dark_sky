"""Step 0: fetch public input data that is reachable without credentials.

- Copernicus DEM GLO-30 (AWS open data), domain + HORIZON_MAX_KM buffer
- ESA WorldCover 2021 v200 10 m (AWS open data), domain bbox (water mask)
- Lorenz World Atlas binary tiles (2022 per brief, 2025 as cross-check)
VIIRS DNB radiance is NOT fetched here: EOG requires a login. See viirs.py.
"""
import gzip
import math
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import rasterio
from rasterio.merge import merge
import requests

from . import config as C


def domain_lv95():
    """Square analysis domain in LV95 (minx, miny, maxx, maxy), snapped to 1 km."""
    from pyproj import Transformer
    x0, y0 = Transformer.from_crs(C.CRS_LL, C.CRS_M, always_xy=True).transform(C.BERN[1], C.BERN[0])
    x0, y0 = round(x0, -3), round(y0, -3)
    h = C.DOMAIN_HALF_KM * 1000
    return x0 - h, y0 - h, x0 + h, y0 + h


def domain_bounds():
    """Lon/lat bounding box of the square domain."""
    from pyproj import Transformer
    minx, miny, maxx, maxy = domain_lv95()
    xs = np.array([minx, maxx, maxx, minx, (minx + maxx) / 2, (minx + maxx) / 2])
    ys = np.array([miny, miny, maxy, maxy, miny, maxy])
    lon, lat = Transformer.from_crs(C.CRS_M, C.CRS_LL, always_xy=True).transform(xs, ys)
    return lon.min(), lat.min(), lon.max(), lat.max()


def dem_bounds():
    minx, miny, maxx, maxy = domain_bounds()
    dlat = C.HORIZON_MAX_KM / 111.0 + 0.02
    dlon = C.HORIZON_MAX_KM / (111.0 * math.cos(math.radians((miny + maxy) / 2))) + 0.02
    return minx - dlon, miny - dlat, maxx + dlon, maxy + dlat


def _cop_url(lat, lon):
    ns = f"N{lat:02d}" if lat >= 0 else f"S{-lat:02d}"
    ew = f"E{lon:03d}" if lon >= 0 else f"W{-lon:03d}"
    name = f"Copernicus_DSM_COG_10_{ns}_00_{ew}_00_DEM"
    return f"https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"


def _download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    with requests.get(url, stream=True, timeout=300) as r:
        if r.status_code == 404:
            return None
        r.raise_for_status()
        tmp = dest.with_suffix(".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
        tmp.rename(dest)
    return dest


def fetch_dem():
    out = C.INTERIM / "dem_wgs84.tif"
    if out.exists():
        return out
    minx, miny, maxx, maxy = dem_bounds()
    tiles = [(la, lo) for la in range(math.floor(miny), math.floor(maxy) + 1)
             for lo in range(math.floor(minx), math.floor(maxx) + 1)]
    d = C.RAW / "copdem"
    d.mkdir(exist_ok=True)
    with ThreadPoolExecutor(6) as ex:
        paths = list(ex.map(lambda t: _download(_cop_url(*t), d / f"cop_{t[0]}_{t[1]}.tif"), tiles))
    paths = [p for p in paths if p]
    srcs = [rasterio.open(p) for p in paths]
    arr, tr = merge(srcs, bounds=(minx, miny, maxx, maxy), nodata=-9999)
    prof = srcs[0].profile | dict(driver="GTiff", height=arr.shape[1], width=arr.shape[2],
                                  transform=tr, compress="deflate", tiled=True, nodata=-9999)
    with rasterio.open(out, "w", **prof) as dst:
        dst.write(arr)
    return out


def fetch_worldcover():
    out = C.INTERIM / "worldcover.tif"
    if out.exists():
        return out
    minx, miny, maxx, maxy = domain_bounds()
    minx, miny, maxx, maxy = minx - 0.02, miny - 0.02, maxx + 0.02, maxy + 0.02
    tiles = {(3 * math.floor(la / 3), 3 * math.floor(lo / 3))
             for la in (miny, maxy) for lo in (minx, maxx)}
    parts = []
    for la, lo in sorted(tiles):
        url = ("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
               f"ESA_WorldCover_10m_2021_v200_N{la:02d}E{lo:03d}_Map.tif")
        parts.append(rasterio.open(url))
    arr, tr = merge(parts, bounds=(minx, miny, maxx, maxy), nodata=0)
    prof = parts[0].profile | dict(driver="GTiff", height=arr.shape[1], width=arr.shape[2],
                                   transform=tr, compress="deflate", tiled=True, nodata=0)
    with rasterio.open(out, "w", **prof) as dst:
        dst.write(arr)
    return out


# ---------------- Lorenz atlas ----------------
def _lorenz_decode(raw):
    """Decode a Lorenz binary tile -> (600,600) brightness ratio (artificial/natural),
    row 0 = southern edge. Mirrors getInfoFromLonLat()/compressed2full() in the
    atlas overlay JS: first value is 2 bytes, then signed 1-byte deltas; column 0
    is a running sum in latitude, each row then a running sum in longitude."""
    a = np.frombuffer(gzip.decompress(raw), dtype=np.int8).astype(np.int64)
    first = 128 * a[0] + a[1]
    d = a[1:].reshape(600, 600).copy()  # d[i, j]; d[0,0] is the low byte of first
    d[0, 0] = 0
    col0 = first + np.cumsum(d[:, 0])
    vals = col0[:, None] + np.concatenate([np.zeros((600, 1), np.int64), np.cumsum(d[:, 1:], axis=1)], axis=1)
    return (5.0 / 195.0) * (np.exp(0.0195 * vals) - 1.0)


def fetch_lorenz(year):
    out = C.INTERIM / f"lorenz{year}_ratio.tif"
    if out.exists():
        return out
    minx, miny, maxx, maxy = dem_bounds()
    tx = sorted({math.floor((lo + 180) / 5) + 1 for lo in (minx, maxx)})
    ty = sorted({math.floor((la + 65) / 5) + 1 for la in (miny, maxy)})
    tx, ty = range(tx[0], tx[-1] + 1), range(ty[0], ty[-1] + 1)
    mosaic = np.zeros((600 * len(ty), 600 * len(tx)), np.float32)
    for jj, y in enumerate(ty):
        for ii, x in enumerate(tx):
            url = f"https://djlorenz.github.io/astronomy/binary_tiles/{year}/binary_tile_{x}_{y}.dat.gz"
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            (C.RAW / f"lorenz_{year}_{x}_{y}.dat.gz").write_bytes(r.content)
            t = _lorenz_decode(r.content)
            r0 = (len(ty) - 1 - jj) * 600
            mosaic[r0:r0 + 600, ii * 600:(ii + 1) * 600] = t[::-1]  # north-up
    # The overlay JS picks index round(120*dx + 1/2) - 1 == floor(120*dx): value k covers
    # [k/120, (k+1)/120) from the tile's SW corner, i.e. area-registered cells.
    west = -180 + 5 * (tx[0] - 1)
    north = -65 + 5 * ty[-1]
    tr = rasterio.transform.from_origin(west, north, 1 / 120, 1 / 120)
    with rasterio.open(out, "w", driver="GTiff", height=mosaic.shape[0], width=mosaic.shape[1], count=1,
                       dtype="float32", crs=C.CRS_LL, transform=tr, compress="deflate") as dst:
        dst.write(mosaic, 1)
    return out


def ratio_to_mpsas(ratio):
    """Lorenz convention: natural zenith 22.0 mag/arcsec^2."""
    return 22.0 - 2.5 * np.log10(1.0 + ratio)


if __name__ == "__main__":
    for fn in (fetch_dem, fetch_worldcover, lambda: fetch_lorenz(2022), lambda: fetch_lorenz(2025)):
        t = time.time()
        print(fn(), f"{time.time() - t:.0f}s", flush=True)
