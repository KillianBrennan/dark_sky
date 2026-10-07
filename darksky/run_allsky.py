"""Driver for the all-sky (> 30 deg) layer with local sources, visibility and building masks.

  .venv/bin/python -m darksky.run_allsky            # central run + layer
  .venv/bin/python -m darksky.run_allsky sens       # fill x0 / x2 sensitivity
  .venv/bin/python -m darksky.run_allsky finalize   # re-apply masks + re-colour, no model run
  .venv/bin/python -m darksky.run_allsky colour     # re-colour only
"""
import json
import sys
import time

import numpy as np
import rasterio
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree

from . import allsky as A
from . import config as C
from . import layers as Ly
from .calibrate import lorenz_at
from .fetch import ratio_to_mpsas
from .lightdome import natural_rel
from .local_sources import build as build_sources
from .terrain import load_dem

AOD_GRID = (0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.25)


def calibrate(H, t250):
    m_obs = ratio_to_mpsas(lorenz_at(t250["x"], t250["y"]))
    wh, wa = np.array([3]), np.array([0])
    rows = []
    for aod in AOD_GRID:
        K = A.kernel(aod)
        B = A._far(np.ascontiguousarray(K), np.log(A.D_KM), A.L.DZ_M, A.L.EPS_DEG, wh, wa, float(C.AZ_BIN),
                   C.SHADOW_LEAK, False, *(np.ascontiguousarray(t250[k], np.float64) for k in ("x", "y", "z")),
                   np.zeros((1, 72, len(t250["beta_d_km"])), np.float32), np.log(t250["beta_d_km"]),
                   H["px"], H["py"], H["pz"], H["pw"], H["ps"], H["cx"], H["cy"], H["cz"], H["cw"], H["cstart"], H["ccount"],
                   H["bx"], H["by"], H["bz"], H["bw"], H["bstart"], H["bcount"],
                   C.SOURCE_RADIUS_KM * 1e3, 15e3, 0.0, C.LOCAL_TAPER_M)[:, 0].astype(np.float64)
        f = lambda lc: np.sqrt(np.mean((ratio_to_mpsas(10 ** lc * B) - m_obs) ** 2))
        r = minimize_scalar(f, bounds=(-5, 25), method="bounded")
        res = ratio_to_mpsas(10 ** r.x * B) - m_obs
        rows.append(dict(aod=aod, log10_scale=float(r.x), rms_mag=float(r.fun), bias_mag=float(res.mean()),
                         p95_abs_mag=float(np.quantile(np.abs(res), 0.95))))
        print(rows[-1], flush=True)
    best = dict(min(rows, key=lambda r: r["rms_mag"]))
    best["grid"] = rows
    best["n_points"] = int(len(m_obs))
    (C.OUT / "calibration_allsky.json").write_text(json.dumps(best, indent=2))
    return best


def viirs_only_sources():
    from .local_sources import _viirs_points
    x, y, w, _, _ = _viirs_points(C.SOURCE_RADIUS_KM * 1000.0)
    return dict(x=x, y=y, w=w, kind=np.zeros(len(x), np.int8)), {"variant": "raw VIIRS pixels as point sources"}


def model(tag, fill_factor, cal=None):
    t = time.time()
    src, st = viirs_only_sources() if tag == "viirs_only" else build_sources(fill_factor)
    H = A.hierarchy(src, load_dem())
    print(f"[{tag}] sources: {len(src['w'])} points, {len(H['cw'])} cells, {len(H['bw'])} blocks ({time.time() - t:.0f}s)", flush=True)
    t250 = np.load(C.INTERIM / "terrain_250m.npz")
    if cal is None:
        cal = calibrate(H, t250)
    K = A.kernel(cal["aod"])
    t = time.time()
    Bf = A.far_term(K, H, (t250["x"], t250["y"], t250["z"]), t250["beta"].astype(np.float32), t250["beta_d_km"])
    print(f"[{tag}] far term: {time.time() - t:.0f}s", flush=True)
    t100 = np.load(C.INTERIM / "terrain_100m.npz")
    t = time.time()
    Bl = A.local_term(K, H, (t100["x"], t100["y"], t100["z"]))
    print(f"[{tag}] local term: {time.time() - t:.0f}s", flush=True)
    np.savez(C.INTERIM / f"allsky_{tag}.npz", far=Bf, local=Bl, aod=cal["aod"], log10_scale=cal["log10_scale"])
    st = {k: v for k, v in st.items()}
    (C.OUT / f"sources_{tag}.json").write_text(json.dumps(st, indent=2, default=float))
    return cal


def building_distance(x, y):
    from .osm import buildings
    b = buildings()
    tree = cKDTree(np.c_[b.x.values, b.y.values])
    d, i = tree.query(np.c_[x, y], k=8, distance_upper_bound=2000.0)
    r = np.where(np.isfinite(d), b.radius_m.values[np.minimum(i, len(b) - 1)], 0.0)
    return np.nanmin(np.where(np.isfinite(d), np.maximum(d - r, 0.0), np.inf), axis=1)


def fields(tag="central"):
    m = np.load(C.INTERIM / f"allsky_{tag}.npz")
    scale = 10 ** float(m["log10_scale"]) * C.LED_FACTOR["V"][0]
    g250 = np.load(C.INTERIM / "terrain_250m.npz")
    g100 = np.load(C.INTERIM / "terrain_100m.npz")
    tr250, sh250 = rasterio.Affine(*g250["transform"]), tuple(g250["shape"])
    tr100 = rasterio.Affine(*g100["transform"])
    far = np.exp(Ly._resample_250_to_100(np.log(np.maximum(m["far"].astype(np.float64), 1e-30)), g250, tr250, sh250, g100, tr100))
    loc = m["local"].astype(np.float64)
    tot = scale * (far + loc)
    wh, _ = A.want_list()
    band = []
    for ih, h in enumerate(A.HS):
        sel = wh == ih
        band.append(natural_rel(h) + tot[:, sel].mean(1))
    mean30 = sum(w * b for w, b in zip(A.BAND_W, band))
    f = {"allsky30_mag": C.NATURAL_ZENITH_MPSAS - 2.5 * np.log10(mean30),
         "zenith_mag": C.NATURAL_ZENITH_MPSAS - 2.5 * np.log10(band[3]),
         "local_share": (scale * loc[:, wh == 3][:, 0]) / tot[:, wh == 3][:, 0]}
    return f, g100, tr100


def tree_fraction(g):
    """ESA WorldCover 2021 tree cover (class 10) share of each 100 m pixel. The Copernicus
    DSM puts the observer on top of the canopy inside large forests, so the horizon test
    alone does not remove them."""
    from rasterio.warp import Resampling, reproject
    with rasterio.open(C.INTERIM / "worldcover.tif") as src:
        trees = (src.read(1) == 10).astype(np.float32)
        s_tr, s_crs = src.transform, src.crs
    sh = tuple(g["shape"])
    out = np.zeros(sh, np.float32)
    reproject(trees, out, src_transform=s_tr, src_crs=s_crs, dst_transform=rasterio.Affine(*g["transform"]),
              dst_crs=C.CRS_M, resampling=Resampling.average)
    return out[g["rows"], g["cols"]]


def run(recalibrate=False):
    cf = C.OUT / "calibration_allsky.json"
    cal = None if recalibrate or not cf.exists() else json.loads(cf.read_text())
    model("central", C.MISSED_FILL_FACTOR, cal)
    finalize()


def finalize():
    """Masks, fields file, layer and summary from the saved central model run."""
    f, g, tr = fields("central")
    sh = tuple(g["shape"])
    vis = (g["horizon5"] < C.VIS_MAX_HORIZON_DEG).mean(1)
    bdist = building_distance(g["x"], g["y"])
    water = Ly._water_fraction(g["x"], g["y"]) > 0.5
    trees = tree_fraction(g)
    keep = ((vis >= C.VIS_MIN_FRACTION) & (bdist >= C.BUILDING_MIN_DIST_M) & (trees < C.MAX_TREE_FRACTION)
            & ~water)
    f["visible_fraction_below10"] = vis
    f["building_dist_m"] = bdist
    f["tree_fraction"] = trees
    f["shown"] = keep.astype(np.float32)
    # previous layer (VIIRS-only sources, 2-point sky sampling) for comparison
    with rasterio.open(C.OUT / "darksky_fields_lv95.tif") as s:
        old = s.read(list(s.descriptions).index("allsky30_mag") + 1)
    f["allsky30_minus_previous_layer"] = f["allsky30_mag"] - old[g["rows"], g["cols"]]
    path = C.OUT / "allsky30_fields_lv95.tif"
    with rasterio.open(path, "w", driver="GTiff", height=sh[0], width=sh[1], count=len(f), dtype="float32",
                       crs=C.CRS_M, transform=tr, nodata=np.nan, compress="deflate", tiled=True) as dst:
        for i, (k, v) in enumerate(f.items(), 1):
            dst.write(Ly._raster(g, v, sh), i)
            dst.set_band_description(i, k)
    shown = np.where(keep, f["allsky30_mag"], np.nan)
    colour_layer(shown, g, tr, sh)
    s = {"shown_fraction_of_land": float(keep.sum() / (~water).sum()),
         "fail_visibility": float(((vis < C.VIS_MIN_FRACTION) & ~water).sum() / (~water).sum()),
         "fail_buildings": float(((bdist < C.BUILDING_MIN_DIST_M) & ~water).sum() / (~water).sum()),
         "fail_forest": float(((trees >= C.MAX_TREE_FRACTION) & ~water).sum() / (~water).sum()),
         "removed_by_forest_only": float(((trees >= C.MAX_TREE_FRACTION) & (vis >= C.VIS_MIN_FRACTION)
                                          & (bdist >= C.BUILDING_MIN_DIST_M) & ~water).sum() / (~water).sum()),
         "allsky30_shown_quantiles": np.nanquantile(shown, [0, .02, .25, .5, .75, .98, 1]).round(2).tolist(),
         "change_vs_previous_layer_shown_quantiles": np.nanquantile(np.where(keep, f["allsky30_minus_previous_layer"], np.nan),
                                                                [.02, .25, .5, .75, .98]).round(3).tolist(),
         "local_share_zenith_shown_median": float(np.nanmedian(np.where(keep, f["local_share"], np.nan)))}
    (C.OUT / "allsky30_summary.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))


BREAKS = [19.75, 20.0, 20.25, 20.5, 20.75, 21.0]   # 7 classes of 0.25 mag


def colour_layer(shown=None, g=None, tr=None, sh=None):
    """Stepped (7-class) RGBA COG + legend; re-colours from the fields file if called alone."""
    if shown is None:
        g = np.load(C.INTERIM / "terrain_100m.npz")
        tr, sh = rasterio.Affine(*g["transform"]), tuple(g["shape"])
        with rasterio.open(C.OUT / "allsky30_fields_lv95.tif") as s:
            dsc = list(s.descriptions)
            mag = s.read(dsc.index("allsky30_mag") + 1)[g["rows"], g["cols"]]
            ok = s.read(dsc.index("shown") + 1)[g["rows"], g["cols"]] > 0
        shown = np.where(ok, mag, np.nan)
    Ly.write_rgba("3_allsky30", shown, g, tr, sh, BREAKS[0], BREAKS[-1], breaks=BREAKS)
    Ly.stepped_legend("3_allsky30", "3  Mean sky brightness above 30°, V (mag/arcsec²), incl. local sources",
                      BREAKS, "Shown where terrain < 10° over ≥ 75 % of the horizon, ≥ 200 m from buildings, < 50 % forest.")


def sensitivity():
    cal = json.loads((C.OUT / "calibration_allsky.json").read_text())
    base, g, _ = fields("central")
    with rasterio.open(C.OUT / "allsky30_fields_lv95.tif") as s:
        shown = s.read(list(s.descriptions).index("shown") + 1)[g["rows"], g["cols"]] > 0
    out = {}
    for tag, ff in (("viirs_only", None), ("fill0", 0.0), ("fill2", 2.0)):
        model(tag, ff, cal)
        f, *_ = fields(tag)
        d = (f["allsky30_mag"] - base["allsky30_mag"])[shown]
        rk = lambda a: np.argsort(np.argsort(a))
        out[tag] = {"delta_mag_quantiles_shown": np.quantile(d, [.02, .5, .98]).round(3).tolist(),
                    "spearman_shown": float(np.corrcoef(rk(base["allsky30_mag"][shown]), rk(f["allsky30_mag"][shown]))[0, 1])}
        print(tag, out[tag], flush=True)
    (C.OUT / "allsky30_sensitivity.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    if sys.argv[1:] == ["sens"]:
        sensitivity()
    elif sys.argv[1:] == ["colour"]:
        colour_layer()
    elif sys.argv[1:] == ["finalize"]:
        finalize()
    else:
        run(recalibrate="--recalibrate" in sys.argv)
