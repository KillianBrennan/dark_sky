"""Step 6: derive the per-target fields on the 100 m grid and write map layers.

Outputs (outputs/):
  darksky_fields_lv95.tif  multi-band float32 data GeoTIFF (EPSG:2056, 100 m), for analysis
  layer_*.tif               RGBA Cloud-Optimised GeoTIFFs (EPSG:2056) for map.geo.admin.ch
                            -> Werkzeuge / Tools > Import > local file
  legend_*.png              colour keys

Fields
  zenith_mag            modelled V zenith brightness, LED-corrected (mag/arcsec^2)
  zenith_mag_noLED      same, Lorenz-equivalent (no LED correction)
  allsky30_mag          mean V brightness over the sky above 30 deg (mag/arcsec^2)
  se25_S_mag, se25_E_mag  V brightness at 25 deg, mean over S (135-225) / E (45-135) sector
  airglow_excess        target 1: artificial NIR (Gen3) skyglow in the darkest usable 30-deg wide
                        window at 20-30 deg elevation within azimuth 45-225 (E..S), in mag
                        relative to the best 1 % of the domain (0 = as dark as the best sites).
                        A 5-deg azimuth bin is usable only if terrain stays below 15 deg there
                        (the bottom of the 15-35 deg airglow observing band).
                        NaN where no usable window exists (terrain blocks E and S).
  airglow_best_az       centre azimuth of that window (deg)
  meteor_rel            target 4: relative visual meteor rate = open-sky fraction above 15 deg
                        x r^(NELM - NELM_pristine), r = 2.5 population index, NELM from
                        allsky30_mag (Schaefer-type relation)
  open15                fraction of the sky above 15 deg not blocked by terrain
  p_above_inversion     P(site above Oct-Feb stratus top), normal(1050 m, 200 m) assumption
  elevation_m           DEM elevation
  lorenz_zenith_mag     Lorenz atlas (reference, independent of terrain)
  blocked_frac_se25_V   fraction of S/E 25-deg artificial V light removed by terrain
"""
import json

import numpy as np
import rasterio
from matplotlib import pyplot as plt
from rasterio.shutil import copy as rio_copy
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy.stats import norm

from . import config as C
from .calibrate import lorenz_at
from .fetch import ratio_to_mpsas  # noqa: F401 (re-exported)
from .lightdome import natural_rel
from .run_model import H_IDX, SE_BINS

# validated sequential blue ramp (dataviz reference palette), light -> dark
BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
ORANGE = "#eb6834"
HATCH = (120, 118, 112)
POP_INDEX = 2.5
WINDOW_BINS = 30 // C.AZ_BIN
MAX_USABLE_HORIZON = 15.0         # terrain must stay below the bottom of the 15-35 deg observing band


def _hex(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)], float)


def ramp_lut(n=256):
    stops = np.array([_hex(h) for h in BLUE])
    pos = np.linspace(0, 1, len(stops))
    t = np.linspace(0, 1, n)
    return np.stack([np.interp(t, pos, stops[:, k]) for k in range(3)], 1).astype(np.uint8)


def nelm(mpsas):
    """Naked-eye limiting magnitude from sky brightness (Schaefer 1990-type fit)."""
    return 7.93 - 5 * np.log10(10 ** (4.316 - mpsas / 5) + 1)


def _grid(npz):
    d = np.load(C.INTERIM / npz)
    tr = rasterio.Affine(*d["transform"])
    return d, tr, tuple(d["shape"])


def _resample_250_to_100(vals250, g250, tr250, shape250, g100, tr100):
    """Bilinear resample of per-pixel 250 m values (n250, k) onto 100 m pixel centres."""
    vals250 = np.atleast_2d(vals250.T).T
    out = np.empty((len(g100["x"]), vals250.shape[1]), np.float32)
    inside = np.zeros(shape250, bool)
    inside[g250["rows"], g250["cols"]] = True
    _, (ri, ci) = distance_transform_edt(~inside, return_indices=True)
    fr = (tr250.f - g100["y"]) / -tr250.e - 0.5
    fc = (g100["x"] - tr250.c) / tr250.a - 0.5
    for k in range(vals250.shape[1]):
        r = np.full(shape250, np.nan, np.float32)
        r[g250["rows"], g250["cols"]] = vals250[:, k]
        r = r[ri, ci]                                   # nearest-fill outside the polygon
        out[:, k] = map_coordinates(r, [fr, fc], order=1, mode="nearest")
    return out


def compute(model_tag="central"):
    cal = json.loads((C.OUT / "calibration.json").read_text())
    scale = 10 ** cal["log10_scale"]
    g250, tr250, sh250 = _grid("terrain_250m.npz")
    g100, tr100, sh100 = _grid("terrain_100m.npz")
    m = np.load(C.INTERIM / f"model_{model_tag}.npz")
    mv = m if "B_V" in m else np.load(C.INTERIM / "model_central.npz")   # NIR-only sensitivity runs
    BV, BN, BbV = mv["B_V"], m["B_NIR"], mv["Bblk_V"]
    led_v, led_n = C.LED_FACTOR["V"][0], C.LED_FACTOR["NIR"][0]
    n_se = len(SE_BINS)

    # assemble 250 m channels: zen, 45-deg ring (8), V25 S/E bins, blocked V25 S/E, NIR mean(20,25,30) S/E bins
    ring = np.arange(0, 360 // C.AZ_BIN, 45 // C.AZ_BIN)
    nir = (BN[:, H_IDX[20], SE_BINS] + BN[:, H_IDX[25], SE_BINS] + BN[:, H_IDX[30], SE_BINS]) / 3
    ch = np.concatenate([BV[:, H_IDX[90], :1], BV[:, H_IDX[45], ring], BV[:, H_IDX[25], SE_BINS],
                         BbV[:, H_IDX[25], SE_BINS], nir], axis=1)
    # interpolate in log space (fields span orders of magnitude)
    ch100 = np.exp(_resample_250_to_100(np.log(np.maximum(ch, 1e-30)), g250, tr250, sh250, g100, tr100))
    zen, r45 = ch100[:, 0], ch100[:, 1:9]
    v25 = ch100[:, 9:9 + n_se]
    v25b = ch100[:, 9 + n_se:9 + 2 * n_se]
    nirse = ch100[:, 9 + 2 * n_se:]

    f = {}
    f["zenith_mag"] = ratio_to_mpsas(scale * led_v * zen)
    f["zenith_mag_noLED"] = ratio_to_mpsas(scale * zen)
    nat45, nat25 = natural_rel(45), natural_rel(25)
    w45 = (np.sin(np.radians(60)) - np.sin(np.radians(30))) / (1 - np.sin(np.radians(30)))
    mean30 = w45 * (nat45 + scale * led_v * r45.mean(1)) + (1 - w45) * (1 + scale * led_v * zen)
    f["allsky30_mag"] = C.NATURAL_ZENITH_MPSAS - 2.5 * np.log10(mean30)
    az_se = SE_BINS * C.AZ_BIN
    S, E = (az_se >= 135) & (az_se <= 225), (az_se >= 45) & (az_se <= 135)
    f["se25_S_mag"] = C.NATURAL_ZENITH_MPSAS - 2.5 * np.log10(nat25 + scale * led_v * v25[:, S].mean(1))
    f["se25_E_mag"] = C.NATURAL_ZENITH_MPSAS - 2.5 * np.log10(nat25 + scale * led_v * v25[:, E].mean(1))
    f["blocked_frac_se25_V"] = v25b.sum(1) / (v25.sum(1) + v25b.sum(1))

    # --- target 1: airglow window ---
    hor = g100["horizon5"][:, SE_BINS]
    usable = hor <= MAX_USABLE_HORIZON
    art = scale * led_n * nirse
    nwin = n_se - WINDOW_BINS + 1
    win_val = np.full((len(art), nwin), np.inf)
    for w in range(nwin):
        sl = slice(w, w + WINDOW_BINS)
        ok = usable[:, sl].all(1)
        win_val[ok, w] = art[ok, sl].mean(1)
    best = win_val.min(1)
    best_w = win_val.argmin(1)
    best = np.where(np.isfinite(best), best, np.nan)
    ref = np.nanquantile(best, 0.01)
    f["airglow_excess"] = 2.5 * np.log10(best / ref)
    f["airglow_best_az"] = np.where(np.isfinite(best), az_se[best_w] + (WINDOW_BINS - 1) * C.AZ_BIN / 2, np.nan)

    # --- target 4: meteors ---
    f["open15"] = g100["open15"]
    f["meteor_rel"] = g100["open15"] * POP_INDEX ** (nelm(f["allsky30_mag"]) - nelm(C.NATURAL_ZENITH_MPSAS))

    # --- inversion ---
    f["elevation_m"] = g100["z"]
    f["p_above_inversion"] = norm.cdf((g100["z"] - C.INVERSION_TOP_MEAN_M) / C.INVERSION_TOP_SD_M)
    f["lorenz_zenith_mag"] = ratio_to_mpsas(lorenz_at(g100["x"], g100["y"]))
    f["drive_min_bin"] = g100["drive_min"]
    # mask lakes/rivers (WorldCover 80, majority within the 100 m pixel): not observable
    water = _water_fraction(g100["x"], g100["y"]) > 0.5
    for k in f:
        if k not in ("elevation_m", "drive_min_bin"):
            f[k] = np.where(water, np.nan, f[k])
    return f, g100, tr100, sh100


def _water_fraction(x, y):
    from pyproj import Transformer
    from scipy.ndimage import uniform_filter
    with rasterio.open(C.INTERIM / "worldcover.tif") as src:
        wc = src.read(1)
        tr = src.transform
    frac = uniform_filter((wc == 80).astype(np.float32), size=9)      # ~ 9 x 10 m
    lon, lat = Transformer.from_crs(C.CRS_M, C.CRS_LL, always_xy=True).transform(x, y)
    r = np.clip(((lat - tr.f) / tr.e).astype(int), 0, wc.shape[0] - 1)
    c = np.clip(((lon - tr.c) / tr.a).astype(int), 0, wc.shape[1] - 1)
    return frac[r, c]


def _raster(g, v, shape, fill=np.nan):
    out = np.full(shape, fill, np.float32)
    out[g["rows"], g["cols"]] = v
    return out


def write_fields(f, g, tr, shape):
    path = C.OUT / "darksky_fields_lv95.tif"
    names = list(f)
    with rasterio.open(path, "w", driver="GTiff", height=shape[0], width=shape[1], count=len(names),
                       dtype="float32", crs=C.CRS_M, transform=tr, nodata=np.nan, compress="deflate",
                       tiled=True) as dst:
        for i, n in enumerate(names, 1):
            dst.write(_raster(g, f[n], shape), i)
            dst.set_band_description(i, n)
    return path


def _hatch(shape, period=6, width=2):
    r, c = np.indices(shape)
    return ((r + c) % period) < width


def stepped_colours(n):
    """n evenly spaced steps of the blue ramp, light -> dark."""
    return [BLUE[i] for i in np.linspace(0, len(BLUE) - 1, n).round().astype(int)]


def write_rgba(name, value, g, tr, shape, vmin, vmax, dark_is_high=True, hatch_nan=None,
               alpha_from=None, breaks=None):
    """Colour `value` on the blue ramp (dark = darkest sky / best) and write an RGBA COG.
    With `breaks` (ascending class boundaries) the ramp is stepped into len(breaks)+1 classes."""
    v = _raster(g, value, shape)
    if breaks is not None:
        cols = np.array([_hex(h) for h in stepped_colours(len(breaks) + 1)], np.uint8)
        k = np.digitize(np.nan_to_num(v, nan=-np.inf), breaks)
        rgb = cols[k if dark_is_high else len(breaks) - k]
    else:
        t = (v - vmin) / (vmax - vmin)
        if not dark_is_high:
            t = 1 - t
        t = np.clip(t, 0, 1)
        lut = ramp_lut()
        idx = np.nan_to_num(t * 255, nan=0).astype(int)
        rgb = lut[idx]
    a = np.where(np.isfinite(v), 255, 0).astype(np.uint8)
    inside = np.zeros(shape, bool)
    inside[g["rows"], g["cols"]] = True
    if hatch_nan is not None:          # boolean per-pixel array: where NaN means "blocked"
        hh = inside & _raster(g, hatch_nan.astype(np.float32), shape, 0).astype(bool) & _hatch(shape)
        rgb[hh] = HATCH
        a[hh] = 230
    if alpha_from is not None:
        rgb[:] = _hex(ORANGE).astype(np.uint8)
        a = (np.clip(np.nan_to_num(_raster(g, alpha_from, shape, 0)), 0, 1) * 220).astype(np.uint8)
        a[~inside] = 0
    tmp = C.INTERIM / f"{name}_tmp.tif"
    with rasterio.open(tmp, "w", driver="GTiff", height=shape[0], width=shape[1], count=4, dtype="uint8",
                       crs=C.CRS_M, transform=tr, photometric="RGB", alpha="UNSPECIFIED") as dst:
        for k in range(3):
            dst.write(rgb[..., k], k + 1)
        dst.write(a, 4)
        dst.colorinterp = [rasterio.enums.ColorInterp.red, rasterio.enums.ColorInterp.green,
                           rasterio.enums.ColorInterp.blue, rasterio.enums.ColorInterp.alpha]
    out = C.OUT / f"layer_{name}.tif"
    rio_copy(tmp, out, driver="COG", compress="DEFLATE", overview_resampling="nearest", blocksize=256)
    tmp.unlink()
    return out


def stepped_legend(name, title, breaks, note=None, fmt="{:.2f}"):
    """Discrete legend: one swatch per class, boundaries labelled between swatches."""
    cols = stepped_colours(len(breaks) + 1)
    fig, ax = plt.subplots(figsize=(4.6, 1.3 if not note else 1.55), dpi=200)
    n = len(cols)
    for i, c in enumerate(cols):
        ax.add_patch(plt.Rectangle((i + 0.03, 0), 0.94, 1, color=c, lw=0))   # 2px-ish surface gap
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks(range(1, n))
    ax.set_xticklabels([fmt.format(b) for b in breaks], fontsize=7, color="#52514e")
    ax.text(0.0, -0.55, "← brighter sky", fontsize=6.5, color="#52514e", ha="left", va="top", transform=ax.transData)
    ax.text(n, -0.55, "darker sky →", fontsize=6.5, color="#52514e", ha="right", va="top", transform=ax.transData)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0, pad=3)
    ax.set_title(title, fontsize=8, loc="left", color="#0b0b0b")
    if note:
        fig.text(0.01, 0.02, note, fontsize=6, color="#52514e")
    fig.tight_layout()
    p = C.OUT / f"legend_{name}.png"
    fig.savefig(p, transparent=False, facecolor="#fcfcfb")
    plt.close(fig)
    return p


def legend(name, title, vmin, vmax, ticks, labels=None, note=None, dark_is_high=True, hatch_label=None):
    lut = ramp_lut() / 255.0
    if not dark_is_high:
        lut = lut[::-1]
    fig, ax = plt.subplots(figsize=(4.6, 1.15 if not note else 1.4), dpi=200)
    ax.imshow(lut[None, :, :], aspect="auto", extent=(vmin, vmax, 0, 1))
    ax.set_yticks([])
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels or [str(t) for t in ticks], fontsize=7, color="#52514e")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0, pad=3)
    ax.set_title(title, fontsize=8, loc="left", color="#0b0b0b")
    if note:
        fig.text(0.01, 0.02, note, fontsize=6, color="#52514e")
    fig.tight_layout()
    p = C.OUT / f"legend_{name}.png"
    fig.savefig(p, transparent=False, facecolor="#fcfcfb")
    plt.close(fig)
    return p
