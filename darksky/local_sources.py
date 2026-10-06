"""Hybrid light-source field: VIIRS downscaled with BFS hectare statistics.

Why: VIIRS DNB has a ~0.5-1 km footprint (plus geolocation scatter), and the masked VNL
product drops pixels below ~0.45 nW/cm^2/sr. Inside Switzerland the BFS hectare grids
(residents STATPOP 2024, jobs STATENT 2023 in FTE, buildings 2024) say *where* light can
be at 100 m. In the domain only ~3 % of residents live where VIIRS sees nothing, but these
farms and hamlets sit precisely in the dark areas that matter for site choice.

Construction (Swiss territory only; abroad VIIRS is used unchanged):
  proxy_h = sum_k c_k v_hk,  v = residents, jobs (FTE), buildings, metres of lit road,
            metres of untagged major road (OSM, see osm.py); c_k >= 0 are the VIIRS
            intensity per unit from a non-negative least-squares fit over 5 km cells
            inside Switzerland, so proxy_h is already in intensity units.
  1. Downscale: each lit VIIRS pixel's intensity is redistributed to hectares within
     600 m with weights proxy_h * exp(-d^2 / 2 sigma^2), sigma = 400 m. If no hectare
     with proxy > 0 lies within 600 m (roads, industry, greenhouses, ski lifts) the light
     stays at the VIIRS pixel centre. Total VIIRS intensity is conserved.
  2. Fill: hectares with proxy > 0 but no lit VIIRS pixel within 600 m get
     MISSED_FILL_FACTOR * proxy_h,
     capped so that the summed fill inside any VIIRS pixel stays below the detection floor
     (otherwise VIIRS would have seen it).
Output: point sources (x, y, z, w) in LV95, w in the same units as viirs.sources().
"""
import json
import math

import numpy as np
import pandas as pd
import rasterio
import requests
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely import contains_xy
from shapely.geometry import shape

from . import config as C
from .fetch import domain_lv95

BOUNDARY_URL = ("https://api3.geo.admin.ch/rest/services/api/MapServer/"
                "ch.swisstopo.swissboundaries3d-land-flaeche.fill/CH?geometryFormat=geojson&sr=2056")
_COLL = {"pop": ("ch.bfs.volkszaehlung-bevoelkerungsstatistik_einwohner", "volkszaehlung-bevoelkerungsstatistik_einwohner_2024"),
         "jobs": ("ch.bfs.betriebszaehlungen-beschaeftigte_vollzeitaequivalente", "betriebszaehlungen-beschaeftigte_vollzeitaequivalente_2023"),
         "bldg": ("ch.bfs.volkszaehlung-gebaeudestatistik_gebaeude", "volkszaehlung-gebaeudestatistik_gebaeude_2024")}


def fetch_bfs():
    d = C.RAW / "bfs"
    d.mkdir(exist_ok=True)
    for k, f in C.BFS_FILES.items():
        p = d / f
        if not p.exists():
            coll, item = _COLL[k]
            r = requests.get(C.BFS_URL.format(coll=coll, item=item, file=f), timeout=300)
            r.raise_for_status()
            p.write_bytes(r.content)
    b = C.RAW / "ch_boundary.json"
    if not b.exists():
        r = requests.get(BOUNDARY_URL, timeout=120)
        r.raise_for_status()
        b.write_text(r.text)
    return d


def swiss_polygon():
    fetch_bfs()
    return shape(json.loads((C.RAW / "ch_boundary.json").read_text())["feature"]["geometry"])


PROXY_VARS = ("pop", "jobs", "bldg", "road_lit", "road_major")


def hectares(margin_m):
    """BFS hectare table + OSM road metres (centre coordinates) within the domain + margin."""
    from .osm import roads_per_hectare
    fetch_bfs()
    minx, miny, maxx, maxy = domain_lv95()
    out = None
    for k, f in C.BFS_FILES.items():
        t = pd.read_csv(C.RAW / "bfs" / f, sep=";", decimal=",", usecols=["E_KOORD", "N_KOORD", "NUMMER"])
        t = t.rename(columns={"NUMMER": k}).groupby(["E_KOORD", "N_KOORD"], as_index=False).sum()
        out = t if out is None else out.merge(t, on=["E_KOORD", "N_KOORD"], how="outer")
    out = out.merge(roads_per_hectare(), on=["E_KOORD", "N_KOORD"], how="outer").fillna(0.0)
    out["x"] = out.E_KOORD + 50.0          # BFS coordinates are the SW corner of the hectare
    out["y"] = out.N_KOORD + 50.0
    m = ((out.x > minx - margin_m) & (out.x < maxx + margin_m) & (out.y > miny - margin_m) & (out.y < maxy + margin_m))
    return out[m].reset_index(drop=True)


def _viirs_points(margin_m):
    """Lit VIIRS pixels as points: LV95 x, y, intensity (same units as viirs.sources), pixel area."""
    from .viirs import clip
    with rasterio.open(clip()) as src:
        rad = src.read(1).astype(np.float64)
        tr = src.transform
    rows, cols = np.nonzero(np.isfinite(rad) & (rad > C.VIIRS_MIN_RADIANCE))
    lon = tr.c + (cols + 0.5) * tr.a
    lat = tr.f + (rows + 0.5) * tr.e
    area = ((abs(tr.a) * math.pi / 180 * C.R_EARTH) * (abs(tr.e) * math.pi / 180 * C.R_EARTH)
            * np.cos(np.radians(lat)) * 1e4)                                  # cm^2
    w = rad[rows, cols] * area * 1e-18
    x, y = Transformer.from_crs(C.CRS_LL, C.CRS_M, always_xy=True).transform(lon, lat)
    return x, y, w, area, tr


def fit_proxy(h, vx, vy, vw, res=5000.0):
    """Non-negative least squares: VIIRS intensity ~ sum_k c_k v_k over res-sized cells
    (Swiss VIIRS pixels and Swiss hectares only). Returns (coef dict, log10 correlation)."""
    from scipy.optimize import nnls
    key = lambda x, y: (np.floor(x / res).astype(np.int64) * 100000 + np.floor(y / res).astype(np.int64))
    kv = pd.Series(vw).groupby(key(vx, vy)).sum().rename("v")
    kh = h.groupby(key(h.x.values, h.y.values))[list(PROXY_VARS)].sum()
    j = kh.join(kv, how="left").fillna(0.0)
    A = j[list(PROXY_VARS)].values
    scale = A.max(0) + 1e-30
    coef, _ = nnls(A / scale, j.v.values)
    coef = coef / scale
    pred = A @ coef
    r = float(np.corrcoef(np.log10(j.v + 1e-6), np.log10(pred + 1e-6))[0, 1])
    return dict(zip(PROXY_VARS, coef)), r


def build(fill_factor=C.MISSED_FILL_FACTOR):
    """Returns dict(x, y, w, kind) of point sources and a stats dict.
    kind: 0 = VIIRS kept in place (abroad, or no proxy nearby), 1 = downscaled to hectare,
    2 = fill for VIIRS-dark inhabited hectare."""
    margin = C.SOURCE_RADIUS_KM * 1000.0
    vx, vy, vw, varea, vtr = _viirs_points(margin)
    ch = swiss_polygon()
    v_in = contains_xy(ch, vx, vy)
    h = hectares(margin_m=margin)
    h = h[contains_xy(ch, h.x.values, h.y.values)].reset_index(drop=True)
    coef, fit_r = fit_proxy(h, vx[v_in], vy[v_in], vw[v_in])
    h["proxy"] = sum(coef[k] * h[k] for k in PROXY_VARS)
    h = h[h.proxy > 0].reset_index(drop=True)
    hx, hy, hp = h.x.values, h.y.values, h.proxy.values

    # 1. downscale Swiss VIIRS pixels onto hectares
    tree = cKDTree(np.c_[hx, hy])
    S = np.zeros(len(h))
    keep = ~v_in
    neigh = tree.query_ball_point(np.c_[vx[v_in], vy[v_in]], r=C.DOWNSCALE_RADIUS_M)
    idx_in = np.nonzero(v_in)[0]
    for i, nb in zip(idx_in, neigh):
        if not nb:
            keep[i] = True
            continue
        nb = np.asarray(nb)
        d2 = (hx[nb] - vx[i]) ** 2 + (hy[nb] - vy[i]) ** 2
        wt = hp[nb] * np.exp(-d2 / (2 * C.DOWNSCALE_SIGMA_M ** 2))
        if wt.sum() <= 0:
            keep[i] = True
            continue
        S[nb] += vw[i] * wt / wt.sum()

    # 2. fill VIIRS-dark inhabited hectares (no lit VIIRS pixel within the radius)
    vtree = cKDTree(np.c_[vx, vy])
    dark = np.array([len(n) == 0 for n in vtree.query_ball_point(np.c_[hx, hy], r=C.DOWNSCALE_RADIUS_M)])
    fill = np.where(dark, fill_factor * hp, 0.0)
    # cap: summed fill within one VIIRS pixel footprint must stay below the detection floor
    lon, lat = Transformer.from_crs(C.CRS_M, C.CRS_LL, always_xy=True).transform(hx, hy)
    pr = np.floor((lat - vtr.f) / vtr.e).astype(np.int64)
    pc = np.floor((lon - vtr.c) / vtr.a).astype(np.int64)
    pix = pr * 1_000_000 + pc
    cap_pix = C.VIIRS_DETECTION_FLOOR * np.median(varea) * 1e-18
    tot = pd.Series(fill).groupby(pix).transform("sum").values
    fill = np.where(tot > cap_pix, fill * cap_pix / np.maximum(tot, 1e-30), fill)

    x = np.concatenate([vx[keep], hx, hx[fill > 0]])
    y = np.concatenate([vy[keep], hy, hy[fill > 0]])
    w = np.concatenate([vw[keep], S, fill[fill > 0]])
    kind = np.concatenate([np.zeros(keep.sum(), np.int8), np.ones(len(hx), np.int8), np.full((fill > 0).sum(), 2, np.int8)])
    sel = w > 0
    stats = dict(viirs_total=float(vw.sum()), viirs_swiss=float(vw[v_in].sum()),
                 kept_in_place_swiss=float(vw[v_in & keep].sum()), downscaled=float(S.sum()),
                 fill_total=float(fill.sum()), n_fill_hectares=int((fill > 0).sum()),
                 n_hectares=int(len(h)), proxy_fit_log_r=fit_r,
                 proxy_share={k: float((coef[k] * h[k]).sum() / h.proxy.sum()) for k in PROXY_VARS},
                 proxy_coef_rel_pop={k: float(coef[k] / coef["pop"]) if coef["pop"] > 0 else None for k in PROXY_VARS})
    return dict(x=x[sel], y=y[sel], w=w[sel], kind=kind[sel]), stats


if __name__ == "__main__":
    import time
    t = time.time()
    src, st = build()
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in st.items()}, f"{time.time() - t:.0f}s")
    print("proxy share:", {k: round(v, 3) for k, v in st["proxy_share"].items()})
    print("coef relative to 1 resident:", {k: (round(v, 4) if v is not None else None) for k, v in st["proxy_coef_rel_pop"].items()})
    print("fill / Swiss VIIRS = %.3f %%" % (100 * st["fill_total"] / st["viirs_swiss"]))
    print("Swiss VIIRS kept in place (no proxy nearby) = %.1f %%" % (100 * st["kept_in_place_swiss"] / st["viirs_swiss"]))
