"""Driver for the all-sky (> 30 deg) layer with local sources, visibility and building masks.

  .venv/bin/python -m darksky.run_allsky            # central run + layer
  .venv/bin/python -m darksky.run_allsky sens       # fill x0 / x2 sensitivity
  .venv/bin/python -m darksky.run_allsky finalize   # re-apply masks + re-colour, no model run
  .venv/bin/python -m darksky.run_allsky colour     # re-colour only
  .venv/bin/python -m darksky.run_allsky glare      # glare term (+ sensitivity), then finalize
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


def glare_term(tag="central", F=C.GARSTANG_F, age=C.CIE_AGE, pig=C.CIE_PIGMENT, fill_factor=C.MISSED_FILL_FACTOR):
    """Veiling luminance from lamps in direct view, averaged over the sky above 30 deg with
    the same solid-angle weights as the sky brightness (natural-zenith units)."""
    from .glare import los_table, veil
    cal = json.loads((C.OUT / "calibration_allsky.json").read_text())
    scale = 10 ** cal["log10_scale"] * C.LED_FACTOR["V"][0]
    src, _ = build_sources(fill_factor)
    H = A.hierarchy(src, load_dem())
    g = np.load(C.INTERIM / "terrain_100m.npz")
    los_f = C.INTERIM / "los_100m.npy"
    if los_f.exists():
        los = np.load(los_f)
    else:
        t = time.time()
        los = los_table(g["x"], g["y"], g["z"])
        np.save(los_f, los)
        print(f"line-of-sight table: {time.time() - t:.0f}s", flush=True)
    t = time.time()
    V, E = veil(H, (g["x"], g["y"], g["z"]), los, scale, cal["aod"], F=F, age=age, pig=pig)
    wh, _ = A.want_list()
    veil30 = sum(w * V[:, wh == ih].mean(1) for ih, w in enumerate(A.BAND_W))
    np.savez(C.INTERIM / f"glare_{tag}.npz", veil30=veil30.astype(np.float32), E=E, F=F, age=age, pig=pig)
    print(f"[{tag}] glare: {time.time() - t:.0f}s", flush=True)
    return veil30


def building_distance(x, y):
    from .osm import buildings
    b = buildings()
    tree = cKDTree(np.c_[b.x.values, b.y.values])
    d, i = tree.query(np.c_[x, y], k=8, distance_upper_bound=2000.0)
    r = np.where(np.isfinite(d), b.radius_m.values[np.minimum(i, len(b) - 1)], 0.0)
    return np.nanmin(np.where(np.isfinite(d), np.maximum(d - r, 0.0), np.inf), axis=1)


def fields(tag="central", glare_tag=None):
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
    gf = C.INTERIM / f"glare_{glare_tag or tag}.npz"
    veil30 = np.load(gf)["veil30"].astype(np.float64) if gf.exists() else np.zeros_like(mean30)
    f = {"allsky30_mag": C.NATURAL_ZENITH_MPSAS - 2.5 * np.log10(mean30 + veil30),
         "allsky30_sky_only_mag": C.NATURAL_ZENITH_MPSAS - 2.5 * np.log10(mean30),
         "glare_veil_ratio": veil30 / mean30,
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
    colour_layer(g=g)
    s = {"shown_fraction_of_land": float(keep.sum() / (~water).sum()),
         "fail_visibility": float(((vis < C.VIS_MIN_FRACTION) & ~water).sum() / (~water).sum()),
         "fail_buildings": float(((bdist < C.BUILDING_MIN_DIST_M) & ~water).sum() / (~water).sum()),
         "fail_forest": float(((trees >= C.MAX_TREE_FRACTION) & ~water).sum() / (~water).sum()),
         "removed_by_forest_only": float(((trees >= C.MAX_TREE_FRACTION) & (vis >= C.VIS_MIN_FRACTION)
                                          & (bdist >= C.BUILDING_MIN_DIST_M) & ~water).sum() / (~water).sum()),
         "allsky30_shown_quantiles": np.nanquantile(shown, [0, .02, .25, .5, .75, .98, 1]).round(2).tolist(),
         "change_vs_previous_layer_shown_quantiles": np.nanquantile(np.where(keep, f["allsky30_minus_previous_layer"], np.nan),
                                                                [.02, .25, .5, .75, .98]).round(3).tolist(),
         "local_share_zenith_shown_median": float(np.nanmedian(np.where(keep, f["local_share"], np.nan))),
         "glare_veil_over_sky_shown_quantiles": np.nanquantile(np.where(keep, f["glare_veil_ratio"], np.nan),
                                                              [.1, .5, .9, .99]).round(3).tolist(),
         "glare_dmag_shown_quantiles": np.nanquantile(np.where(keep, f["allsky30_mag"] - f["allsky30_sky_only_mag"], np.nan),
                                                     [.01, .1, .5]).round(3).tolist(),
         "breaks": breaks_from_unmasked(f["allsky30_mag"])}
    (C.OUT / "allsky30_summary.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))


def breaks_from_unmasked(mag, n=7, step=0.01):
    """Logarithmic (halving) class breaks over the unmasked field: from the brightest end the
    classes hold 50 %, 25 %, 12.5 %, ... of the pixels and the darkest class the remainder,
    so colour resolution is concentrated in the dark tail. Rounded to `step` mag."""
    q = np.nanquantile(mag, 1.0 - 0.5 ** np.arange(1, n))
    return [round(float(np.round(v / step) * step), 2) for v in q]


def colour_layer(shown=None, g=None, tr=None, sh=None):
    """Stepped 7-class RGBA COGs + legends from the fields file: masked layer_3 and fully
    unmasked layer_3b. Breaks are halving quantiles of the *unmasked* field."""
    g = np.load(C.INTERIM / "terrain_100m.npz") if g is None else g
    tr, sh = rasterio.Affine(*g["transform"]), tuple(g["shape"])
    with rasterio.open(C.OUT / "allsky30_fields_lv95.tif") as s:
        dsc = list(s.descriptions)
        mag = s.read(dsc.index("allsky30_mag") + 1)[g["rows"], g["cols"]]
        ok = s.read(dsc.index("shown") + 1)[g["rows"], g["cols"]] > 0
    br = breaks_from_unmasked(mag)
    title = "3  Sky above 30°, V (mag/arcsec²), incl. local sources and glare"
    Ly.write_rgba("3_allsky30", np.where(ok, mag, np.nan), g, tr, sh, br[0], br[-1], breaks=br)
    Ly.stepped_legend("3_allsky30", title, br,
                      "Shown where terrain < 10° over ≥ 75 % of the horizon, ≥ 200 m from buildings, < 50 % forest.")
    Ly.write_rgba("3b_allsky30_unmasked", mag, g, tr, sh, br[0], br[-1], breaks=br)
    Ly.stepped_legend("3b_allsky30_unmasked", title, br,
                      "No masks. Classes hold 50, 25, 12.5, 6.25, 3.1, 1.6, 1.6 % of the square (brightest to darkest).")
    return br


def glare_sensitivity():
    """Glare under other assumptions: horizontal emission F = 0.30, observer age 25 and 60."""
    base, g, _ = fields("central")
    with rasterio.open(C.OUT / "allsky30_fields_lv95.tif") as s:
        shown = s.read(list(s.descriptions).index("shown") + 1)[g["rows"], g["cols"]] > 0
    out = {}
    rk = lambda a: np.argsort(np.argsort(a))
    for tag, kw in (("glare_F030", dict(F=0.30)), ("glare_age25", dict(age=25.0)), ("glare_age60", dict(age=60.0))):
        glare_term(tag, **kw)
        f, *_ = fields("central", glare_tag=tag)
        d = (f["allsky30_mag"] - base["allsky30_mag"])[shown]
        out[tag] = {"delta_mag_quantiles_shown": np.quantile(d, [.02, .5, .98]).round(3).tolist(),
                    "spearman_shown": float(np.corrcoef(rk(base["allsky30_mag"][shown]), rk(f["allsky30_mag"][shown]))[0, 1])}
        print(tag, out[tag], flush=True)
    (C.OUT / "allsky30_glare_sensitivity.json").write_text(json.dumps(out, indent=2))


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
    elif sys.argv[1:] == ["glare"]:
        glare_term()
        finalize()
        glare_sensitivity()
    else:
        run(recalibrate="--recalibrate" in sys.argv)
