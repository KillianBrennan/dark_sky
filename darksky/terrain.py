"""DEM handling + step 3: horizon profiles by ray-casting.

The Copernicus GLO-30 DSM is reprojected once to LV95 (EPSG:2056) at DEM_RES_M.
LV95 scale error is <0.1 % within ~150 km of Bern, negligible for horizon angles.

Horizon: for each site and azimuth, march outward with step max(DEM_RES, s/150)
to HORIZON_MAX_KM, bilinear-sample terrain, and take the max apparent elevation
    alt = atan((z - z0 - s^2 / (2 k R)) / s),   k = 4/3 (standard refraction).
Observer eye height 1.6 m. Copernicus is a DSM, so forest/buildings are included
in the profile, which is realistic for near-field obstruction.

Also returns, per site and 5-degree azimuth sector, the terrain *shadow angle*
seen from a source at distance D (needed by the light-dome blocking model):
    beta(D) = max_{s<D} atan((z(s) - z(D) - curvature) / (D - s))
i.e. how high above its own horizontal a ground source at distance D must emit
for its light to clear the terrain between it and the observer.
"""
import math

import numpy as np
import rasterio
from numba import njit, prange
from rasterio.warp import Resampling, calculate_default_transform, reproject

from . import config as C

BETA_D_KM = np.geomspace(1.0, 300.0, 28)  # source distances for shadow-angle table
BETA_AZ_STEP = 5


def dem_lv95():
    out = C.INTERIM / "dem_lv95.tif"
    if out.exists():
        return out
    src_p = C.INTERIM / "dem_wgs84.tif"
    with rasterio.open(src_p) as src:
        tr, w, h = calculate_default_transform(src.crs, C.CRS_M, src.width, src.height, *src.bounds,
                                               resolution=C.DEM_RES_M)
        prof = src.profile | dict(crs=C.CRS_M, transform=tr, width=w, height=h, dtype="float32",
                                  nodata=-9999, compress="deflate", tiled=True, BIGTIFF="IF_SAFER")
        with rasterio.open(out, "w", **prof) as dst:
            reproject(rasterio.band(src, 1), rasterio.band(dst, 1), resampling=Resampling.bilinear,
                      dst_nodata=-9999)
    return out


_DEM = None


def load_dem():
    """3x3 median filter removes single-pixel DSM spikes (telecom masts, pylons such as the
    Chasseral tower) that otherwise appear as 20-degree walls, while preserving ridgelines."""
    global _DEM
    if _DEM is None:
        from scipy.ndimage import median_filter
        with rasterio.open(dem_lv95()) as src:
            z = src.read(1).astype(np.float32)
            tr = src.transform
        z[z < -1000] = np.nan
        _DEM = (median_filter(z, size=3), tr)
    return _DEM


@njit(cache=True, fastmath=True)
def _bilinear(z, fx, fy):
    ny, nx = z.shape
    if fx < 0 or fy < 0 or fx >= nx - 1 or fy >= ny - 1:
        return np.nan
    i = int(fy)
    j = int(fx)
    dy = fy - i
    dx = fx - j
    return ((z[i, j] * (1 - dx) + z[i, j + 1] * dx) * (1 - dy)
            + (z[i + 1, j] * (1 - dx) + z[i + 1, j + 1] * dx) * dy)


@njit(parallel=True, cache=True, fastmath=True)
def _raycast(z, x0, y0, res, xs, ys, z0s, n_az, max_m, re_eff, beta_d, beta_az_step, do_beta):
    n = xs.shape[0]
    hor = np.full((n, n_az), -90.0, dtype=np.float32)
    hor_dist = np.zeros((n, n_az), dtype=np.float32)
    nb_az = 360 // beta_az_step
    nd = beta_d.shape[0]
    beta = np.full((n if do_beta else 1, nb_az, nd), -90.0, dtype=np.float32)
    # pre-compute step distances
    s = res
    count = 0
    while s < max_m:
        count += 1
        s += max(res, s / 150.0)
    dists = np.empty(count)
    s = res
    for k in range(count):
        dists[k] = s
        s += max(res, s / 150.0)
    for p in prange(n):
        zo = z0s[p]
        prof = np.empty(count)
        for a in range(n_az):
            az = math.radians(a * 360.0 / n_az)
            sa = math.sin(az)
            ca = math.cos(az)
            best = -1e9
            bestd = 0.0
            for k in range(count):
                d = dists[k]
                fx = (xs[p] + d * sa - x0) / res - 0.5
                fy = (y0 - (ys[p] + d * ca)) / res - 0.5
                zz = _bilinear(z, fx, fy)
                if np.isnan(zz):
                    prof[k] = np.nan
                    continue
                # height relative to observer's tangent plane (refraction via k-factor)
                u = zz - d * d / (2.0 * re_eff)
                prof[k] = u
                t = (u - zo) / d
                if t > best:
                    best = t
                    bestd = d
            hor[p, a] = math.degrees(math.atan(best))
            hor_dist[p, a] = bestd
            # shadow angle table, only on sector-centre azimuths
            if do_beta and (a * 360 // n_az) % beta_az_step == 0:
                ib = (a * 360 // n_az) // beta_az_step
                for j in range(nd):
                    D = beta_d[j]
                    # source ground height = terrain at distance D (index of nearest sample)
                    kd = 0
                    while kd < count - 1 and dists[kd] < D:
                        kd += 1
                    us = prof[kd]
                    if np.isnan(us):
                        us = 400.0 - D * D / (2.0 * re_eff)  # off-DEM: assume lowland source
                    b = -1e9
                    for k in range(kd):
                        if not np.isnan(prof[k]):
                            t = (prof[k] - us) / (D - dists[k])
                            if t > b:
                                b = t
                    # angle measured at the source relative to *its* local horizontal:
                    # tilt by D / R (source vertical leans away from observer)
                    beta[p, ib, j] = math.degrees(math.atan(b)) - math.degrees(D / re_eff)
    return hor, hor_dist, beta


def horizons(xs, ys, z0, az_step=C.AZ_STEP_DEG, with_beta=True):
    """xs, ys in LV95 metres; z0 ground elevation. Returns horizon[n, 360/az_step] (deg),
    horizon distance (m), shadow angle table beta[n, 72, len(BETA_D_KM)] (deg; a dummy
    (1, 72, nD) array if with_beta is False)."""
    z, tr = load_dem()
    re_eff = C.K_REFRACTION * C.R_EARTH
    return _raycast(z, tr.c, tr.f, tr.a, np.asarray(xs, float), np.asarray(ys, float),
                    np.asarray(z0, float) + 1.6, int(360 / az_step), C.HORIZON_MAX_KM * 1000.0,
                    re_eff, BETA_D_KM * 1000.0, BETA_AZ_STEP, with_beta)


def open_fraction(hor, min_el=15.0):
    """Fraction of the solid angle above `min_el` not blocked by terrain."""
    h = np.maximum(hor, min_el)
    return ((1 - np.sin(np.radians(np.minimum(h, 90)))).mean(axis=1)) / (1 - math.sin(math.radians(min_el)))
