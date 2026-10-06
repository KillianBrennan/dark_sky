"""Step 5: calibrate the zenith model against the Lorenz atlas.

Terrain blocking is switched OFF here because the Lorenz/Cinzano-type atlas does not
model mountain screening; comparing like with like. Free parameters: overall scale C
(absorbs radiance->flux units and kernel normalisation) and aerosol optical depth
(controls the distance falloff shape). Fit in magnitude space:
    m_model = 22 - 2.5 log10(1 + C * B_zenith(AOD))
NOT identifiable from this comparison (and therefore set as priors in config): the LED
correction (Lorenz is VIIRS-derived and shares the blindness), the Garstang F term and
the terrain-blocking split (Lorenz has no terrain).
"""
import json

import numpy as np
import rasterio
from pyproj import Transformer
from scipy.optimize import minimize_scalar

from . import config as C
from . import lightdome as L
from .fetch import fetch_lorenz, ratio_to_mpsas
from .viirs import sources

AOD_GRID = (0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35)


def lorenz_at(x, y, year=C.LORENZ_YEAR):
    lon, lat = Transformer.from_crs(C.CRS_M, C.CRS_LL, always_xy=True).transform(x, y)
    with rasterio.open(fetch_lorenz(year)) as src:
        return np.array([v[0] for v in src.sample(zip(lon, lat))])


def fit_scale(B, m_obs):
    def rms(logc):
        return np.sqrt(np.mean((ratio_to_mpsas(10 ** logc * B) - m_obs) ** 2))
    r = minimize_scalar(rms, bounds=(-5, 25), method="bounded")
    return r.x, r.fun


def run():
    d = np.load(C.INTERIM / "terrain_250m.npz")
    near, far = sources()
    ratio = lorenz_at(d["x"], d["y"])
    m_obs = ratio_to_mpsas(ratio)
    want = np.ones((1, 1), bool)
    rows = []
    for aod in AOD_GRID:
        K = L.build_kernel(aod=aod, zenith_only=True, bands=("V",))
        B, _ = L.site_brightness(K[0], (d["x"], d["y"], d["z"]), d["beta"][:1].astype(np.float32),
                                 d["beta_d_km"], near, far, want, use_terrain=False)
        logc, rms = fit_scale(B[:, 0, 0], m_obs)
        res = ratio_to_mpsas(10 ** logc * B[:, 0, 0]) - m_obs
        rows.append(dict(aod=aod, log10_scale=logc, rms_mag=rms, bias_mag=float(res.mean()),
                         r=float(np.corrcoef(np.log10(10 ** logc * B[:, 0, 0]), np.log10(ratio))[0, 1])))
        print(rows[-1], flush=True)
    rows = [{k: float(v) for k, v in r.items()} for r in rows]
    best = dict(min(rows, key=lambda r: r["rms_mag"]))
    # residual diagnostics at the best fit
    K = L.build_kernel(aod=best["aod"], zenith_only=True, bands=("V",))
    B, _ = L.site_brightness(K[0], (d["x"], d["y"], d["z"]), d["beta"][:1].astype(np.float32),
                             d["beta_d_km"], near, far, want, use_terrain=False)
    m_mod = ratio_to_mpsas(10 ** best["log10_scale"] * B[:, 0, 0])
    res = m_mod - m_obs
    bins = [(18, 20), (20, 21), (21, 21.5), (21.5, 22.1)]
    best["residual_by_lorenz_mag"] = {f"{a}-{b}": dict(n=int(((m_obs >= a) & (m_obs < b)).sum()),
                                                        mean=float(res[(m_obs >= a) & (m_obs < b)].mean()),
                                                        sd=float(res[(m_obs >= a) & (m_obs < b)].std()))
                                      for a, b in bins if ((m_obs >= a) & (m_obs < b)).sum() > 0}
    best["n_points"] = int(len(m_obs))
    best["p95_abs_residual_mag"] = float(np.quantile(np.abs(res), 0.95))
    best["lorenz_year"] = C.LORENZ_YEAR
    best["grid"] = rows
    (C.OUT / "calibration.json").write_text(json.dumps(best, indent=2))
    np.save(C.INTERIM / "calib_residual_250m.npy", res)
    return best


if __name__ == "__main__":
    b = run()
    print(json.dumps({k: v for k, v in b.items() if k != "grid"}, indent=2))
